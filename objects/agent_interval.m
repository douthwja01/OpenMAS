%% INTERVAL ANALYSIS TOOL SET (agent_interval.m) %%%%%%%%%%%%%%%%%%%%%%%%%
% This agent is designed to contain all the basic interval methods used in
% interval themed agents that are generic.

% Author: James A. Douthwaite & Allen De Frietas

classdef agent_interval < agent
    %  INITIALISE THE PROPERTIES OF THE FORMATION CONTROLLER
    properties

    end
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods
        % Constructor
        function obj = agent_interval(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent(varargin);                                           % Get the supercalss
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    
    % /////////////////////// SENSOR INTERACTIONS /////////////////////////
    methods  
        % SENSOR MODEL - CAMERA & RANGE FINDER INTERVALS (2D & 3D)
        function [rangeBox,azimuthBox,elevationBox,alphaBox] = GetSensorMeasurementIntervals(obj,obstacleData) 
            % This function takes the simulation data and calculates the
            % spherical position and radius that would otherwise be sensed 
            % by the system.
            % INPUTS:
            % obstacleData - The object observation structure
            % OUTPUTS:
            % range        - The interval in the obstacles range
            % azimuth      - The interval in the angular position in the XY plane
            % elevation    - The interval in the angular position in the XZ plane
            % alpha        - The interval in the obstacles angular width in the azimuth

            % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
            [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.GetSensorMeasurements(obstacleData);
            % DEFINE THE MEASUREMENT INTERVALS
            rangeBox     = midrad(measuredRange,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_rangeFinder);
            if inf(rangeBox) < 0
                rangeBox = infsup(inf(rangeBox),sup(rangeBox));            % Omit the negative range readings
            end
            azimuthBox   = midrad(measuredAzimuth,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_camera);
            elevationBox = midrad(measuredElevation,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_camera);
            alphaBox     = midrad(measuredAlpha,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_camera);
        end
        % SENSOR MODEL - LOCAL IMU & GPS INTERVALS (2D & 3D)
        function [p_i,v_i,r_i] = GetAgentMeasurementIntervals(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            % INPUTS:
            % .SENSORS
            % - sigma_position     - The position measurement standard deviation
            % - sigma_velocity     - The velocity measurement standard deviation
            % - confidenceInterval - The number of standard deviations (typically 3)
            % localState           - The current state of the object 
            % .VIRTUAL.radius      - The radius of the object
            
            % Input sanity check
            assert(isstruct(obj.SENSORS) && ~isempty(obj.SENSORS),'The SENSOR structure is absent.');
            assert(isnumeric(obj.localState), 'The state of the object is invalid.');
            assert(~any(isnan(obj.localState)),'The state contains NaNs');
            
            % Depending on the problem dimension
            if obj.VIRTUAL.is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices = 1:2;
                velocityIndices = 4:5;
            end
            % Get the position and velocity + perturb the state
            positionAUncertainty = obj.SENSORS.sigma_position*randn(numel(positionIndices),1);  % Zero estimate of position
            velocityAUncertainty = obj.SENSORS.sigma_velocity*randn(numel(velocityIndices),1);  % Zero estimate of velocity
            % Generate the intervals centered around that perturbed value
            p_i = obj.localState(positionIndices,1) + midrad(positionAUncertainty,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_position);
            v_i = obj.localState(velocityIndices,1) + midrad(velocityAUncertainty,obj.SENSORS.confidenceAssumption*obj.SENSORS.sigma_velocity); 
            % Assume the knowledge of the radius is absolute
            r_i = obj.VIRTUAL.radius;
        end
        
        % GET OBJECT DATA STRUCTURE (INTERVAL SENSORS)
        function [objectData1] = GetObjectUpdateFromSensors(obj,dt,observedObject)
            % This function overrides the 'agent' function that employs the
            % .SENSOR parameters to emulate the data recieved from an
            % on-board sensor system.
            
            % 1. Get prior knowledge of the object
            % Get previous state
            objectData0 = obj.GetLastEntryFromMemory(observedObject.objectID);
            
            % Create new memory structure
            firstEntry = 0;
            if numel(objectData0) == 0 
                fprintf('Is first record');
                firstEntry = 1;
                % Load template structure for regularity
                [objectData1] = obj.GetMemoryEntryTemplate();
                % META parameters (static)
                objectData1.name     = observedObject.name;
                objectData1.objectID = observedObject.objectID;
                objectData1.type     = observedObject.type;            
            else
                objectData1 = objectData0;
            end
            % META parameters (Dynamic)
            objectData1.priority = observedObject.priority;
            objectData1.TTC      = 1/observedObject.timeToCollision;
            
            % 2. Get new knowledge of the object
            % Get the current sensor readings for the entity
            [d_measured,psi_measured,theta_measured,alpha_measured] = obj.GetSensorMeasurements(observedObject);
            
            % Process the measurements -> state estimation?
            
            % Store the measurements 
            objectData1.range     = d_measured;
            objectData1.azimuth   = psi_measured;
            objectData1.elevation = theta_measured;
            objectData1.alpha     = alpha_measured;
            
                
            
            % Equivalent position
            [measuredPosition] = obj.GetCartesianFromSpherical(d_measured,psi_measured,theta_measured);
            % Equivalent radius
            [measuredRadius]   = obj.GetObjectRadius(d_measured,alpha_measured);
            
            objectData1.position = measuredPosition;
            objectData1.radius   = measuredRadius;
            
            if firstEntry
                return
            end
            
            measuredVelocity = (objectData1.position - objectData0.position)/dt;
            objectData1.velocity = measuredVelocity;
            

        end
    end
    % //////////////// INTERVAL SENSOR PARAMETER IMPORTS  /////////////////
    methods (Static)
        % GET EMPTY MEMORY STRUCTURE
        function [entry] = GetMemoryEntryTemplate()
            % Create emptry memory structure
            entry = struct(...
                'name',[],...
                'objectID',0,...
                'type',[],...
                'position',[],...
                'velocity',[],...
                'radius',[],...
                'range',[],...
                'azimuth',[],...
                'elevation',[],...
                'alpha',[],...
                'TTC',[],...
                'priority',[]);
        end
        % SENSOR MODEL - REPRESENTATIVE SENSING
        function [SENSORS] = GetCustomSensorParameters()
            % This function is designed to populate the SENSOR structure
            % with representative sensing parameters for the assembly of
            % the sensing intervals.
            % BUILD THE SENSOR FACTOR
            SENSORS = struct(...
                'range',inf,...             % Assume the agent has perfect environmental knowledge (m)
                'sigma_position',0.1,...    % Accurate to within 0.5m
                'sigma_velocity',0.1,...    % Accurate to within 0.1m/s
                'sigma_rangeFinder',0.1,... % Accurate to within 0.1m
                'sigma_camera',5.208E-5,... % One pixel in a 1080p image
                'sampleFrequency',inf,...   % Object has perfect precision
                'confidenceAssumption',3);  % Zero standard deviations
        end
        % SENSOR MODEL - PERFECT SENSING
        function [SENSORS] = GetDefaultSensorParameters()
            % This function is designed to populate the SENSOR field with
            % perfect sensor capabilities.
            % BUILD THE SENSOR FACTOR
            SENSORS = struct(...
                'range',inf,...             % Assume the agent has perfect environmental knowledge (m)
                'sigma_position',0.0,...    % Perfect position measurement
                'sigma_velocity',0.0,...    % Perfect velocity measurement
                'sigma_rangeFinder',0.0,... % Perfect range acquistion
                'sigma_camera',0.0,...      % Infinte resolution
                'sampleFrequency',inf,...   % Infinite sampling rate
                'confidenceAssumption',0);  % The interval width is zero
        end 
    end
    % /////////////////////// INTERVAL MATH TOOLS /////////////////////////
    methods (Static, Access = public)
        % GEt LINEAR STATE ESTIMATION
        function [p1,v1,a1] = linearStateEstimation(dt,p0,v0,p1)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
                      
           dp = (p1 - p0); % Position \delta
           
           v1 = dp/dt;     % Defines the average velocity

           a1 = (v1 - v0)/dt;
           
%            p1 = p0 + v1
%            a1 = (v1 - v0)/dt
           
           % NO PREVIOUS VELOCITY RECORDED
           if any(isnan(v0))
              return
           end 
           
           

        end
        % INTERVAL HEADING
        function [v_unit] = heading(v)
            
            % Input sanity check
            if ~isintval(v)
                v_unit = v/norm(v); % If the variable is not an interval
                return
            end

            v_sup = sup(v);
            v_inf = inf(v); 
            v_sup = v_sup/norm(v_sup);
            v_inf = v_inf/norm(v_inf);
            
            v_unit = infsup(-1,1);
            % Get the heading interval
            for i = 1:size(v,1)
                if v_inf(i) > v_sup(i)
                    v_unit(i,1) = infsup(v_sup(i),v_inf(i));
                else
                    v_unit(i,1) = infsup(v_inf(i),v_sup(i));
                end
            end
        end
        % INTERVAL UNIT VECTOR
        function [v_int] = unit(v,zeroOmitted)
            % INPUT HANDLING
            if ~exist('zeroOmitted','var')
                zeroOmitted = 0;
            end
            % CHECK IF NOT INTERVAL TYPE
            if ~isintval(v)
                % If the input value is not an interval
                v_int = v/norm(v);
                return
            end
            % IF ZERO NOT IN RANGE
            if zeroOmitted
                for i = 1:size(v)
                    s = 0; % Collect the square terms
                    for j = 1:size(v)
                        if i == j
                            continue
                        end
                        s = s + v(j)^2;
                    end
                    x = v(i);
                    
                    % COMPUTE THE INTERVAL NORM WHERE 0 IS NOT INCLUDED
                    n(i) = sign(inf(x))/sqrt(1+(s/x^2));
                end
            else
                % ZERO IS PART OF THE RANGE ( GENERALISED )
                for i = 1:size(v,1)
                    s = 0; % Collect the square terms
                    for j = 1:size(v,1)
                        if i == j
                            continue
                        end
                        s = s + sqr(v(j));
                    end
                    x = v(i);
                    
                    % DEFINE THE UPPER/LOWER INTERSECTIONS
                    upperbound = infsup(0,inf);
                    lowerbound = infsup(-inf,0);
                    xi_pos = intersect(x,upperbound);
                    xi_neg = intersect(x,lowerbound);               
                    
                    % DEFINE THE BOUNDS OF THE INTERVAL                                                          
                    ni_upper =  1 / (sqrt(1 + (s/(xi_pos).^2)));
                    ni_lower = -1 / (sqrt(1 + (s/(xi_neg).^2)));
                    
                    if isnan(ni_upper)
                        n(i) = ni_lower;
                    elseif isnan(ni_lower)
                        n(i) = ni_upper;
                    else
                        n(i) = hull(ni_upper,ni_lower);
                    end
                    % IF NO UNION IS FOUND, THEN IT CONTAINS A ZERO
                    if isnan(ni_lower) && isnan(ni_upper)
                        n(i) = infsup(-1,1); % The catcha clause for zero
                    end
                end
            end
            v_int = n;
            % CONFIRM FORMATTING OF THE VECTOR
            v_int = reshape(v_int,size(v));
        end
        % INTERVAL DOT PRODUCT
        function [v_int] = dot(v_a,v_b)
            v_int = v_a(1)*v_b(1) + v_a(2)*v_b(2) + v_a(3)*v_b(3);
        end
        % INTERVAL CROSS PRODUCT
        function [v_int] = cross(v_a,v_b)
            v_int = [v_a(2)*v_b(3)- v_a(3).*v_b(2); 
                     v_a(3)*v_b(1)- v_a(1).*v_b(3); 
                     v_a(1)*v_b(2)- v_a(2).*v_b(1)]; 
        end
        % INTERVAL NORM VECTOR
        function [n_int] = norm(v)
            % This function calculates the norm of an interval vector
            % CHECK IF NOT INTERVAL TYPE
            if ~isintval(v)
                % If the input value is not an interval
                n_int = norm(v);
                return
            end
%             sqrVector = intVector.^2;
%             sumVector = sum(sqrVector);
%             normInt = sqrt(sumVector);
            
%             offset = 0.2E-2;
%             if inf(normInt) == 0
%                 normInt = infsup(offset ,sup(normInt));
%             end

            infVal = inf(v).^2;
            supVal = sup(v).^2;
            sqrMatrix = zeros(size(infVal,1),2);
            for ind = 1:length(infVal)
                if infVal(ind) >= supVal(ind)
                    sqrMatrix(ind,:) = [supVal(ind),infVal(ind)];
                else
                    sqrMatrix(ind,:) = [infVal(ind),supVal(ind)];
                end
            end
            columnSum = sum(sqrMatrix);
            sqrtMatrix = sqrt(columnSum);
            n_int = infsup(sqrtMatrix(1),sqrtMatrix(2));     
        end
    end
    % //////////////////////// DRAWING METHODS ////////////////////////////
    methods (Static)
        % DRAW A VECTOR GEOMETRIC INTERVAL
        function [patchObj] = drawIntervalVector(u,v,plotLogical,colour)
            % Input arguements
            if nargin < 4
                colour = 'b';
            end
            if nargin < 3
                plotLogical = logical(false);
            end
            % Input sanity check
            assert(size(u,2) == 1 && size(v,2) == 1,'The provided vectors must be column vectors.');
            assert(islogical(plotLogical) == 1,'The plot indicator must be a logical.');
            assert(ischar(colour) || isnumeric(colour),'The colour specifier must be a string or RBG vector.');
                    
            
            % Dimensions of the interval vectors
            pointDensity = 2;
            
            %%%% INTERVAL A %%%%
            if isintval(u)
                dimA = size(u,1);
                minimumsA = inf(u);
                maximumsA = sup(u);            
                % CONSTANTS
                numPoints = pointDensity^dimA;
                % CONTAINERS
                dimMultipliers = zeros(dimA,1);
                dimensionalDistribution = zeros(dimA,pointDensity);
                for dim = 1:dimA
                    dimensionalDistribution(dim,:) = [minimumsA(dim),maximumsA(dim)];  % Matrix with the dimensional distribution
                    dimMultipliers(dim) = numPoints/pointDensity^dim;                                       % Dimension tessilation indicies
                end     
                dimMultipliers = fliplr(dimMultipliers);

                % GENERATE THE N-DIMENSIONAL POINT LIST
                intervalExtentsA = zeros(dimA,size(dimensionalDistribution,2)^dimA);
                for dim = 1:dimA
                    insertVector = [];
                    for i = 1:pointDensity
                       distributionIter = repmat(dimensionalDistribution(dim,i),1,dimMultipliers(dim));
                       insertVector = horzcat(insertVector,distributionIter);
                    end
                    no_copies = numPoints/size(insertVector,2);
                    intervalExtentsA(dim,:) = repmat(insertVector,1,no_copies);
                end 
            else
                % Just append the point directly
                intervalExtentsA = u;
            end
            % we assume u is the origin, v is the velocity 
            w = v + u;  
            
            %%%% INTERVAL B %%%%
            if isintval(w)
                dimB = size(w,1);
                minimumsB = inf(w);
                maximumsB = sup(w);
                % CONSTANTS
                numPoints = pointDensity^dimB;
                % CONTAINERS
                dimMultipliers = zeros(dimB,1);
                dimensionalDistribution = zeros(dimB,pointDensity);
                for dim = 1:dimB
                    dimensionalDistribution(dim,:) = [minimumsB(dim),maximumsB(dim)];  % Matrix with the dimensional distribution
                    dimMultipliers(dim) = numPoints/pointDensity^dim;                                       % Dimension tessilation indicies
                end     
                dimMultipliers = fliplr(dimMultipliers);

                % GENERATE THE N-DIMENSIONAL POINT LIST
                intervalExtentsB = zeros(dimB,size(dimensionalDistribution,2)^dimB);
                for dim = 1:dimB
                    insertVector = [];
                    for i = 1:pointDensity
                       distributionIter = repmat(dimensionalDistribution(dim,i),1,dimMultipliers(dim));
                       insertVector = horzcat(insertVector,distributionIter);
                    end
                    no_copies = numPoints/size(insertVector,2);
                    intervalExtentsB(dim,:) = repmat(insertVector,1,no_copies);
                end 
            else
                % Just append the point directly
                intervalExtentsB = w;
            end            
            % Assemble the volume
            K = convhull(...
                [intervalExtentsA(1,:),intervalExtentsB(1,:)]',...
                [intervalExtentsA(2,:),intervalExtentsB(2,:)]',...
                [intervalExtentsA(3,:),intervalExtentsB(3,:)]',...
                'simplify',logical(true));
            
            % Assemble a sudo-patch object
            patchObj = struct('faces',K,...
                'vertices',[intervalExtentsA,intervalExtentsB]');
            
            % Plot if requested
            if plotLogical
                % Centroid vector
                q = quiver3(mid(u(1)),mid(u(2)),mid(u(3)),...
                            mid(v(1)),mid(v(2)),mid(v(3)),'k');
                q.AutoScaleFactor = 1;
                q.LineWidth = 2;
                % 
                patch(gca,...
                    'Faces',patchObj.faces,...
                    'Vertices',patchObj.vertices,...
                    'FaceColor',colour,...
                    'FaceAlpha',0.2);
            end
        end
        % DRAW A GEOMETRIC INTERVAL
        function [patchObj] = drawInterval(u,plotLogical,colour)
            % Input arguements
            if nargin < 3
                colour = 'b';
            end
            if nargin < 2
                plotLogical = logical(false);
            end
            % Input sanity check
            assert(size(u,2) == 1,'The provided vectors must be column vectors.');
            assert(islogical(plotLogical) == 1,'The plot indicator must be a logical.');
            assert(ischar(colour) || isnumeric(colour),'The colour specifier must be a string or RBG vector.');
            
            % Dimensions of the interval vectors
            pointDensity = 2;
            %%%% INTERVAL A %%%%
            if isintval(u)
                dimA = size(u,1);
                minimumsA = inf(u);
                maximumsA = sup(u);            
                % CONSTANTS
                numPoints = pointDensity^dimA;
                % CONTAINERS
                dimMultipliers = zeros(dimA,1);
                dimensionalDistribution = zeros(dimA,pointDensity);
                for dim = 1:dimA
                    dimensionalDistribution(dim,:) = [minimumsA(dim),maximumsA(dim)];  % Matrix with the dimensional distribution
                    dimMultipliers(dim) = numPoints/pointDensity^dim;                                       % Dimension tessilation indicies
                end     
                dimMultipliers = fliplr(dimMultipliers);

                % GENERATE THE N-DIMENSIONAL POINT LIST
                intervalExtents = zeros(dimA,size(dimensionalDistribution,2)^dimA);
                for dim = 1:dimA
                    insertVector = [];
                    for i = 1:pointDensity
                       distributionIter = repmat(dimensionalDistribution(dim,i),1,dimMultipliers(dim));
                       insertVector = horzcat(insertVector,distributionIter);
                    end
                    no_copies = numPoints/size(insertVector,2);
                    intervalExtents(dim,:) = repmat(insertVector,1,no_copies);
                end 
            else
                % Just append the point directly
                intervalExtents = u;
            end
            % Assemble the volume
            K = convhull(...
                intervalExtents(1,:)',...
                intervalExtents(2,:)',...
                intervalExtents(3,:)',...
                'simplify',logical(true));
            % Assemble a sudo-patch object
            patchObj = struct(...
                'faces',K,...
                'vertices',intervalExtents');
            
            % Plot if requested
            if plotLogical
                patch(gca,...
                    'Faces',patchObj.faces,...
                    'Vertices',patchObj.vertices,...
                    'FaceColor',colour,...
                    'FaceAlpha',0.2);
            end
        end
    end
    % ///////////// BACKGROUND FUNCTIONS (INTERVAL OVERRIDES) /////////////
    methods    
        % TARGET HEADING (INTERVAL)
        function [headingVector]   = GetTargetHeading(obj,targetObject)
            
            % Input check
            if nargin > 1
                headingVector = obj.heading(targetObject.position);               
            elseif ~isempty(obj.targetWaypoint)
                headingVector = obj.heading(obj.targetWaypoint.position);
            else
                % Default heading
                if obj.VIRTUAL.is3D
                    headingVector = [1;0;0];
                else
                    headingVector = [1;0];
                end
            end
        end
        % WAY-POINT CONDITION (INTERVAL))
        function [waypointLogical] = GetTargetCondition(obj)
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            waypointLogical = 0 > norm(mid(obj.targetWaypoint.position)) - (mid(obj.targetWaypoint.radius) + obj.VIRTUAL.radius);
        end
        % MEMORY SORTER (INTERVAL)
        function [reorderedMemory] = SortMemoryByField(obj,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            if ~any(strcmp(fieldnames(obj.memory),field))
                error('Field does not belong to the memory structure.')                
            end
            
            % Handle interval parameters
            if isintval([obj.memory.(field)])
                % Reorder the memory structure based on fieldname
                [~,ind] = sort(mid([obj.memory.(field)]),2,'descend');     % Ordered indices of the object IDs
            else
                [~,ind] = sort([obj.memory.(field)],2,'descend');          % Ordered indices of the object IDs
            end
            % Sort the memory structure 
            reorderedMemory = obj.memory(ind);
        end 
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]