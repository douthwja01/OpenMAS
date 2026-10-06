%% INTERVAL ANALYSIS TOOL SET (agent_interval.m) %%%%%%%%%%%%%%%%%%%%%%%%%
% This agent is designed to contain all the basic interval methods used in
% interval themed agents that are generic.

% Author: James A. Douthwaite & Allen De Frietas

classdef agent_interval < agent
    %  INITIALISE THE PROPERTIES OF THE FORMATION CONTROLLER
    properties
        finiteInf = 1E-11;
    end
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods
        % Constructor
        function this = agent_interval(varargin)
            
            % Get the superclass
            this@agent(varargin);
            
            % Ensure the toolbox is loaded
            IntLab();   % Attempt to launch intlab
            
            % Assign defaults
            this = this.SetBufferSize(3);  
            this.DYNAMICS = this.CreateDYNAMICS();
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Main (default for interval options)
        function this = main(this,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % this      - The updated project
            
            % PLOT AGENT FIGURE
            plotFlag = 1;
            plottingAgent = 1;
            plottedAgent = 2;
            if this.objectID == plottingAgent && plotFlag == 1
                % Plot the measurements in coordinate space
                figure1 = figure(1);
                ax1 = gca;
                hold on; grid on;
                axis equal;
                title(ax1,'Coordinate space');
                xlabel(ax1,'x_{m}'); ylabel(ax1,'y_{m}'); zlabel(ax1,'z_{m}');
                
                % Plot the measurements in temporal space
                figure2 = figure(2);
                ax2 = gca;
                hold on; grid on;
                axis equal;
                title(ax2,'Temporal space');
                xlabel(ax2,'t_{s}'); ylabel(ax2,'x_{-}');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});
            
            % Only plot if we are plotting agent
            if this.objectID == plottingAgent  && plotFlag == 1
                % Check input data is sensed correctly
                %this.GetMeasurementPlot(plottedAgent);
                
                [t] = this.GetTrajectoryByObjectID(plottedAgent,'time');
                [y] = this.GetTrajectoryByObjectID(plottedAgent,'position');
                
                [fh] = this.plotSamples(ax1,y(1:3,:));
                [fh] = this.plotSamplesByDimension(ax2,t,y);
            end
            
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            [headingVector] = this.GetTargetHeading();
            desiredVelocity = mid(headingVector*this.v_nominal);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % [figureHandle] = this.GetObjectScene(gcf);
            
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [this] = this.Controller(ENV.dt,desiredVelocity);
        end
    end
   
    %% /////////////// INTERVAL SENSOR PARAMETER IMPORTS  /////////////////
    methods (Static)
        % SENSOR MODEL - REPRESENTATIVE SENSING
        function [SENSORS] = GetCustomSensorParameters()
            % This function is designed to populate the SENSOR structure
            % with representative sensing parameters for the assembly of
            % the sensing intervals.
            
            % BUILD THE SENSOR STRUCTURE
            SENSORS = struct();
            SENSORS.range = 10;                % Assume the agent has perfect environmental knowledge (m)
            SENSORS.sigma_position = 0.5;      % Accurate to within 0.5m
            SENSORS.sigma_velocity = 0.1;      % Accurate to within 0.1m/s
            SENSORS.sigma_rangeFinder = 0.1;   % Accurate to within 0.1m
            SENSORS.sigma_camera = 5.208E-5;   % One pixel in a 1080p image
            SENSORS.sampleFrequency = inf;     % Object has perfect precision
            SENSORS.confidenceAssumption = 3;  % Zero standard deviations
        end
        % SENSOR MODEL - PERFECT SENSING
        function [SENSORS] = GetDefaultSensorParameters()
            % This function is designed to populate the SENSOR field with
            % perfect sensor capabilities.
            
            % BUILD THE SENSOR STRUCTURE
            SENSORS = struct();
            SENSORS.range = inf;               % Assume the agent has perfect environmental knowledge (m)
            SENSORS.sigma_position = 0.0;      % Perfect position measurement
            SENSORS.sigma_velocity = 0.0;      % Perfect velocity measurement
            SENSORS.sigma_rangeFinder = 0.0;   % Perfect range acquistion
            SENSORS.sigma_camera = 0;          % Infinite resolution
            SENSORS.sampleFrequency = inf;     % Infinite sampling rate
            SENSORS.confidenceAssumption = 0;  % The interval width is zero
        end
    end
    %% ////////////////////// SENSOR INTERACTIONS /////////////////////////
    % The sensor methods employed in the interval agents are specific to
    % methods in this class. The GetAgentUpdate function is used as a
    % comparative 'ideal' case.
    methods
        % Update from sensors defined in the 'SENSORS' (USING SPHERICAL DATA & SENSOR MODEL)
        function [observedObject] = SensorModel(this,dt,observedObject)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observationSet  - The full observed object structure
            % OUTPUTS:
            % this             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            % Input sanity check
            if isempty(observedObject)
                return
            end
            
            % ///// SENSOR MODEL - CAMERA & RANGE INTERVALS (2D & 3D) /////
            % Get measurements from a camera (Noisy or ideal depending on
            % the SENSOR definition.
            [psi_j,theta_j,alpha_j] = this.GetCameraMeasurements(observedObject);
            % Convert the camera measurements to interval boxes
            azimuthBox   = midrad(psi_j,  this.SENSORS.confidenceAssumption*this.SENSORS.sigma_camera);
            elevationBox = midrad(theta_j,this.SENSORS.confidenceAssumption*this.SENSORS.sigma_camera);
            alphaBox     = midrad(alpha_j,this.SENSORS.confidenceAssumption*this.SENSORS.sigma_camera);
            
            % Get the range estimate
            d_j = this.GetRangeFinderMeasurements(observedObject);
            
            % Convert the range-finder measurement to an interval box
            rangeBox = midrad(d_j,this.SENSORS.confidenceAssumption*this.SENSORS.sigma_rangeFinder);
            % Range cannot be less than zero
            rangeBox = intersect(infsup(0,inf),rangeBox);
            
            % /////////////// CALCULATE RADIUS FROM WIDTH /////////////////
            radiusBox = this.GetRadiusFromAngularWidth(rangeBox,alphaBox);
            
            % /////////////////////// STORE DATA //////////////////////////
            % Override the observed object with the sensors intervals
            observedObject.range     = rangeBox;
            observedObject.heading   = azimuthBox;
            observedObject.elevation = elevationBox;
            observedObject.width     = alphaBox;
            observedObject.radius    = radiusBox;
            
            % ///////////// CALCULATE THE CARTESIAN POSITION //////////////
            positionBox = GetCartesianFromSpherical(rangeBox,azimuthBox,elevationBox);
            if ~this.Is3D
                positionBox = positionBox(1:2,1); % Must work for 2D and 3D
            end
            % /////////////////////////////////////////////////////////////
            
            % Up to this point, the measurement intervals are stored in
            % 'observedObject'. From here the position and velocity
            % intervals for the object must be estimated.
                    
            % Define the position and velocity
            observedObject.position = positionBox;
            observedObject.velocity = observedObject.velocity + midrad(0,3*this.SENSORS.sigma_velocity);
        end
        % SENSOR MODEL - LOCAL IMU & GPS INTERVALS (2D & 3D)
        function [p_i,v_i,r_i] = GetAgentMeasurements(this)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in this.localState).
            % INPUTS:
            % .SENSORS
            % - sigma_position     - The position measurement standard deviation
            % - sigma_velocity     - The velocity measurement standard deviation
            % - confidenceInterval - The number of standard deviations (typically 3)
            % localState           - The current state of the object
            % radius               - The radius of the object
            
            % Input sanity check
            assert(isstruct(this.SENSORS) && ~isempty(this.SENSORS),'The SENSOR structure is absent.');
            assert(isnumeric(this.localState),'The state of the object is invalid.');
            assert(~any(isnan(this.localState)),'The state contains NaNs');
            
            % Depending on the problem dimension
            if this.Is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices = 1:2;
                velocityIndices = 4:5;
            end
            % Get the position and velocity + perturb the state
            p_uncertainty = this.SENSORS.sigma_position*randn(numel(positionIndices),1);  % Zero estimate of position
            v_uncertainty = this.SENSORS.sigma_velocity*randn(numel(velocityIndices),1);  % Zero estimate of velocity
            
            % Generate the intervals centered around that perturbed value
            p_i = this.localState(positionIndices,1) + midrad(p_uncertainty,this.SENSORS.confidenceAssumption*this.SENSORS.sigma_position);
            v_i = this.localState(velocityIndices,1) + midrad(v_uncertainty,this.SENSORS.confidenceAssumption*this.SENSORS.sigma_velocity);
            
            % Assume the knowledge of the radius is absolute
            r_i = this.radius;
        end
    end
    
    %% /////////////////////// STATE ESTIMATIION //////////////////////////
    methods (Static)
        % ESTIMATE POSITION
        function [p] = linKin_position(dt,p0,v0)
            p = p0 + v0*dt;
        end
        % ESIMATE DIFFERENTIAL
        function [v] = linKin_velocity(dt,p0,p1)
            v = (p1-p0)/dt;
        end
    end
    
    %% /////////////////////// DYNAMICS & CONTROL ////////////////////////

    % /////////////////////// INTERVAL MATH TOOLS /////////////////////////
    methods (Static, Access = public)
        % Euler heading from vector heading (2D & 3D) - INTERVAL OVERRIDE
        function [lambda,theta] = GetVectorHeadingAngles(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS:
            % V - The current unity velocity vector
            % U - The unit correction vector
            % OUTPUTS:
            % lambda - The azimuth angle (2D)
            % theta  - The pitch angle   (3D)
            
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];                                        % Reject the vertical elements
            % GET THE LINE OF SIGHT ANGLE
            rotationAxis = cross(Vh,Uh);
            
            % HANDLE 3D CASE (ELEVATION & LOS)
            if numel(V) == 3 && numel(U) == 3
                % GET THE LINE OF SIGHT ANGLE
                lambda = sign(mid(rotationAxis(3)))*acos(dot(Vh,Uh)/norm(Vh));  % Get the angle, signed by the direction of its cross product
                % GET THE ELEVATION ANGLE
                theta = atan2(U(3),norm(Uh));
            else
            % HANDLE 2D CASE (LOS)
                % GET THE LINE OF SIGHT ANGLE
                lambda = sign(mid(rotationAxis(3)))*acos(dot(Vh,Uh)/norm(Vh));
                % GET THE SUDO ELEVATION ANGLE
                theta = 0;
            end
        end
        % Interval heading
        function [v_unit] = iheading(v)
            
            % Input sanity check
            if ~isintval(v)
                v_unit = v/norm(v); % If the variable is not an interval
                return
            end
            
            v_sup = sup(v);
            v_inf = inf(v);
            v_sup_unit = v_sup/norm(v_sup);  % The vector to the top left corner
            v_inf_unit = v_inf/norm(v_inf);  % The vector to the bottom right corner
            
            v_unit = infsup(-1,1);
            % Get the heading interval
            for i = 1:size(v,1)
                if v_inf_unit(i) > v_sup_unit(i)
                    v_unit(i,1) = infsup(v_sup_unit(i),v_inf_unit(i));
                else
                    v_unit(i,1) = infsup(v_inf_unit(i),v_sup_unit(i));
                end
            end
        end
        % Interval unit vector
        function [v_int]  = iunit(v)
            
            % Calculate the interval norm
            vnorm = agent_interval.inorm(v);
            if inf(vnorm) == 0
                vnorm = infsup(1E-12,sup(vnorm));
            end
            % Divide the interval vector by the interval norm
            v_int = agent_interval.idivide(v,vnorm);  
        end
        % Interval dot product
        function [v_int]  = idot(v_a,v_b)
            v_int = 0;
            for i = 1:size(v_a,1)
                v_int = v_int + v_a(i)*v_b(i);
            end
        end
        % Interval determinant
        function [v_det]  = idet(v_a,v_b)
            v_det = v_a(1)*v_b(2) - v_b(1)*v_a(2);
        end
        % Interval cross product
        function [v_int]  = icross(v_a,v_b)
            v_int = [...
                v_a(2)*v_b(3)- v_a(3).*v_b(2);
                v_a(3)*v_b(1)- v_a(1).*v_b(3);
                v_a(1)*v_b(2)- v_a(2).*v_b(1)];
        end
        % Interval vector norm
        function [n_int]  = inorm(v)
            % This function calculates the norm of an interval vector
            % CHECK IF NOT INTERVAL TYPE
            
            % Sanity check #1
            if ~isintval(v)         % If the input value is not an interval
                n_int = norm(v);
                return
            end
            % Sanity check #2
            if sum(iszero(v)) == length(v)
               n_int = intval(1E-8); 
               return
            end
            
            infVal = inf(v).^2;
            supVal = sup(v).^2;
            sqrMatrix = zeros(size(infVal,1),2);
            for ind = 1:length(infVal)
                if infVal(ind) > supVal(ind)
                    sqrMatrix(ind,:) = [supVal(ind),infVal(ind)];
                else
                    sqrMatrix(ind,:) = [infVal(ind),supVal(ind)];
                end
            end
            columnSum = sum(sqrMatrix,1);
            sqrtMatrix = sqrt(columnSum);
            
            % Catch zero
            if sqrtMatrix(1) == 0
                n_int = infsup(1E-8,sqrtMatrix(2));
            else
                n_int = infsup(sqrtMatrix(1),sqrtMatrix(2));
            end
        end
        % Interval division (zero catch)
        function [v_int]  = idivide(v_a,v_b)
            % Divide interval (with zero handling)
            %function [z] = divide_patch(y,x)
            % this function is a patch added to intlab for division operation.
            %in intlab, a division x/y and where y contains 0 can lead to NAN values
            % for instance
            %                 xx= infsup(0,2), yy= infsup(0,4), zz=xx/yy
            %                 intval xx =
            %                 [   0.00000000000000,   2.00000000000000]
            %                 intval yy =
            %                 [   0.00000000000000,   4.00000000000000]
            %                 intval zz =
            %                                 NaN
            
            % this is a correct spirit for intlab but this leads to incorrect values
            % for contraction algorithms.
            
            v_int = v_a./v_b;
            ii = find(inf(v_b).*sup(v_b) <= 0);
            infinity = 999999999999;
            for i=ii
                if inf(v_b(i)) < 0 && sup(v_b(i)) > 0
                    % x(i) contains 0'
                    v_int(i) = infsup(-infinity, infinity);
                elseif inf(v_b(i)) == 0 && sup(v_b(i)) ~= 0
                    % If 
                    dummy = infsup(1/sup(v_b(i)), infinity);
                    v_int(i) = v_a(i)*dummy;
                elseif inf(v_b(i)) ~=0 && sup(v_b(i))==0
                    dummy = infsup(-infinity, inf(v_b(i)));
                    v_int(i) = v_a(i)*dummy;
                end
            end
        end
    end
    % //////////////////////// DRAWING METHODS ////////////////////////////
    % These methods are for the plotting of interval measurements/samples.
    methods
        % Plot the incommming sensor readings
        function [ax] = plotMeasurementsByID(this,plotID)
            
            if nargin < 2
                plotID = this.objectID + 1;
            end
            
            % Check input data is sensed correctly
            if ~this.isInMemory(plotID)
                warning('Object not plotted, object %d not known to object %d',plotID,this.objectID);
            end
            
            ax = gca;
            % The time vector
            t_j = this.GetTrajectoryByObjectID(plotID,'time');
            % The range vector
            r_j = this.GetTrajectoryByObjectID(plotID,'range');
            plot(ax,t_j',sup(r_j)','b*','LineStyle','-');
            plot(ax,t_j',inf(r_j)','b*','LineStyle','-.');
            
            % The heading vector
            h_j = this.GetTrajectoryByObjectID(plotID,'heading');
            plot(ax,t_j',sup(h_j)','ro','LineStyle','-');
            plot(ax,t_j',inf(h_j)','ro','LineStyle','-.');
            
            % The elevation vector
            e_j = this.GetTrajectoryByObjectID(plotID,'elevation');
            plot(ax,t_j',sup(e_j)','go','LineStyle','-');
            plot(ax,t_j',inf(e_j)','go','LineStyle','-.');
            
            % The angular width vector
            a_j = this.GetTrajectoryByObjectID(plotID,'width');
            plot(ax,t_j',sup(a_j)','mo','LineStyle','-');
            plot(ax,t_j',inf(a_j)','mo','LineStyle','-.');
        end
    end
    methods (Static)
        % PLOT INTERVAL SERIES - TOGETHER
        function [fh] = plotSamples(ax,y)
            %  This function plots a series of intervals. It is assumed that
            %  x and y are interval vectors
            
            if nargin < 2
                y = ax;
                ax = gca;
            end
            
            [c,d] = size(y);
            
            % Reorient data to regularise plotting
            if c < d
                % Same number of rows
                maxDim = d;
                y = y';
            else % c > d
                % Same number of columns
                maxDim = c;
            end
            % Assume the vectors are concatinated along the longest dimension
            yInf = y;
            ySup = y;
            if isintval(y)
                yInf = NaN(d,c);
                ySup = NaN(d,c);
                for i = 1:maxDim
                    yInf(i,:) = inf(y(i,:));
                    ySup(i,:) = sup(y(i,:));
                end
            end
            
            % With the interval vertices extracted, begin drawing lines
            % between them
            hold on;
            % For each vector
            intDim = size(yInf,2);
            switch intDim
                case 3
                    % Across the complete timeseries (most resent is first
                    for i = 1:maxDim
                        % Check a valid volume exists
                        if any(isnan(yInf(i,:))) || any(isnan(ySup(i,:)))
                            continue
                        end
                        
                        % Define the points in coordinate space
                        vertices(1,:) = [ySup(i,1),ySup(i,2),ySup(i,3)];
                        vertices(2,:) = [ySup(i,1),ySup(i,2),yInf(i,3)];
                        vertices(3,:) = [ySup(i,1),yInf(i,2),yInf(i,3)];
                        vertices(4,:) = [ySup(i,1),yInf(i,2),ySup(i,3)];
                        vertices(5,:) = [yInf(i,1),yInf(i,2),yInf(i,3)];
                        vertices(6,:) = [yInf(i,1),yInf(i,2),ySup(i,3)];
                        vertices(7,:) = [yInf(i,1),ySup(i,2),ySup(i,3)];
                        vertices(8,:) = [yInf(i,1),ySup(i,2),yInf(i,3)];
                        
                        if size(unique(vertices,'rows'),1) > 1
                            K = convhull(vertices(:,1),vertices(:,2),vertices(:,3),'simplify',logical(true));
                            
                            % Assemble a sudo-patch object
                            patchObj = struct('faces',K,'vertices',vertices);
                            % Draw the lines
                            patch(ax,...
                                'Faces',patchthis.faces,'Vertices',patchthis.vertices,...
                                'FaceColor','b','FaceAlpha',0.05);
                        end
                        %                         try
                        %                             K = convhull(vertices(:,1),vertices(:,2),vertices(:,3),'simplify',logical(true));
                        %                         catch
                        %                             warning('Interval volume plot error.');
                        %                             continue
                        %                         end
                        
                        %                         % Assemble a sudo-patch object
                        %                         patchObj = struct('faces',K,'vertices',vertices);
                        %                         % Draw the lines
                        %                         patch(ax,...
                        %                             'Faces',patchthis.faces,'Vertices',patchthis.vertices,...
                        %                             'FaceColor','b','FaceAlpha',0.05);
                        
                        
                        % Check a valid volume exists
                        if i == maxDim || any(isnan(yInf(i+1,:))) || any(isnan(ySup(i+1,:)))
                            continue
                        end
                        
                        % The interval data
                        inf1 = yInf(i,:);
                        sup1 = ySup(i,:);
                        inf2 = yInf(i+1,:);
                        sup2 = ySup(i+1,:);
                        
                        for k = 1:numel(inf1)
                            if inf1(k) < inf2(k)
                                tempMins(1,k) = inf1(k);
                            else
                                tempMins(1,k) = inf2(k);
                            end
                            if sup1(k) > sup2(k)
                                tempMaxs(1,k) = sup1(k);
                            else
                                tempMaxs(1,k) = sup2(k);
                            end
                        end
                        vertices(1,:) = [tempMaxs(1,1),tempMaxs(1,2),tempMaxs(1,3)];
                        vertices(2,:) = [tempMaxs(1,1),tempMaxs(1,2),tempMins(1,3)];
                        vertices(3,:) = [tempMaxs(1,1),tempMins(1,2),tempMins(1,3)];
                        vertices(4,:) = [tempMaxs(1,1),tempMins(1,2),tempMaxs(1,3)];
                        vertices(5,:) = [tempMins(1,1),tempMins(1,2),tempMins(1,3)];
                        vertices(6,:) = [tempMins(1,1),tempMins(1,2),tempMaxs(1,3)];
                        vertices(7,:) = [tempMins(1,1),tempMaxs(1,2),tempMaxs(1,3)];
                        vertices(8,:) = [tempMins(1,1),tempMaxs(1,2),tempMins(1,3)];
                        
                        K = convhull(vertices(:,1),vertices(:,2),vertices(:,3),'simplify',logical(true));
                        % Assemble a sudo-patch object
                        patchObj = struct('faces',K,'vertices',vertices);
                        % Draw the li
                        patch(ax,'Faces',patchthis.faces,'Vertices',patchthis.vertices,...
                            'FaceColor','k','FaceAlpha',0.00,'LineStyle','-.');
                        
                    end
                case 2
                    % Draw the sequencial squares
                    for i = 1:maxDim
                        % Define the points in coordinate space
                        x1 = [yInf(i,1);yInf(i,2)];
                        x2 = [ySup(i,1);yInf(i,2)];
                        x3 = [ySup(i,1);ySup(i,2)];
                        x4 = [yInf(i,1);ySup(i,2)];
                        % Get the square coordinates
                        coords = GetSquare(x1,x2,x3,x4);
                        % Plot the assemble lines
                        fh(1) = plot(ax,coords(:,1),coords(:,2));
                        set(fh(1),'LineStyle','-','Color',[0,0,1]);
                        % Shade the tube
                        fh(2) = fill3(ax,coords(:,1),coords(:,2),coords(:,3),'b','FaceAlpha',0.05);
                    end
                otherwise
                    warning('Accepts only 2D/3D point vectors');
            end
            fh = 0;
        end
        % PLOT INTERVAL SERIES- SEPARATE DIMENSIONS
        function [fh] = plotSamplesByDimension(ax,t,y)
            %  This function plots a series of intervals. It is assumed that
            %  x and y are interval vectors
            
            if nargin < 3
                y = t;
                t = ax;
                ax = gca;
            end
            
            [a,b] = size(t);
            [c,d] = size(y);
            
            % Sanity check
            assert( a == c || b == d,"There must at least one dimension of comparable dimension.");
            assert( ~isintval(t),"Expecting the numeric time vector.");
            
            % Reorient data to regularise plotting
            if a == c
                % Same number of rows
                tmax = a;
            else
                % Same number of columns
                tmax = b;
                t = t'; y = y';
            end
            
            yInf = y;
            ySup = y;
            if isintval(y)
                yInf = NaN(d,c);
                ySup = NaN(d,c);
                for i = 1:tmax
                    yInf(i,:) = inf(y(i,:));
                    ySup(i,:) = sup(y(i,:));
                end
            end
            
            % With the interval vertices extracted, begin drawing lines
            % between them
            hold on;
            % Across the complete timeseries (most resent is first
            for j = 1:tmax
                if isnan(t(j))
                    continue  % Likely the sample is yet to be taken
                end
                % For each vector
                for i = 1:size(yInf,2)
                    % There is no historical
                    if j == tmax
                        continue
                    end
                    
                    % The interval data
                    inf1 = yInf(j,i);   sup1 = ySup(j,i);
                    inf2 = yInf(j+1,i); sup2 = ySup(j+1,i);
                    
                    % Draw enclosing tube
                    x1 = [t(j);inf1]; x2 = [t(j+1);inf2];
                    x3 = [t(j+1);sup2]; x4 = [t(j);sup1];
                    coords = GetSquare(x1,x2,x3,x4);
                    
                    % Plot the assemble lines
                    fh(1) = plot3(ax,coords(:,1),coords(:,2),coords(:,3));
                    set(fh(1),'LineStyle','-','Color',[0,0,1]);
                    
                    % Shade the tube
                    fh(2) = fill3(ax,coords(:,1),coords(:,2),coords(:,3),'b','FaceAlpha',0.05);
                    
                    % Determine the inferior inf
                    if inf1 < inf2
                        x1 = [t(j);inf1];
                        x2 = [t(j+1);inf1];
                    else
                        x1 = [t(j);inf2];
                        x2 = [t(j+1);inf2];
                    end
                    % Determine the superior sup
                    if sup1 > sup2
                        x3 = [t(j+1);sup1];
                        x4 = [t(j);sup1];
                    else
                        x3 = [t(j+1);sup2];
                        x4 = [t(j);sup2];
                    end
                    coords = GetSquare(x1,x2,x3,x4);    % Square coords
                    % Draw the encapsulating square
                    fh(3) = plot3(ax,coords(:,1),coords(:,2),coords(:,3));
                    set(fh(3),'LineStyle','-.','Color','k','LineWidth',0.1);
                    % Draw the interface between the regions
                    fh(4) = line(ax,[t(j);t(j)],[inf1;sup1],'Color','b','LineWidth',1.0);
                end
            end
            hold off;
        end
    end
    % ///////////// BACKGROUND FUNCTIONS (INTERVAL OVERRIDES) /////////////
    methods
        % TARGET HEADING (INTERVAL)
        function [headingVector] = GetTargetHeading(this,targetObject)
            % This is an interval override for the target heading
            % resolution. 
            
            % Input check
            if nargin > 1
                targetPosition = this.GetLastMeasurement(targetObject.objectID,'position');
            elseif ~isempty(this.targetWaypoint)
                targetPosition = this.targetWaypoint.position(:,this.targetWaypoint.sampleNum);
            else
                % Default heading
                if this.Is3D()
                    targetPosition = [1;0;0];
                else
                    targetPosition = [1;0];
                end
            end
            % Get the (interval hardened) vector heading of the target
            headingVector = this.iheading(targetPosition);
        end
        % WAY-POINT CONDITION (INTERVAL))
        function [targetLogical] = GetTargetCondition(this)
            % Get the current measurements
            currentPosition = this.targetWaypoint.position(:,this.targetWaypoint.sampleNum);
            currentRadius   = this.targetWaypoint.radius(:,this.targetWaypoint.sampleNum);
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            targetLogical = 0 > this.inorm(mid(currentPosition)) - (mid(currentRadius) + this.radius);
        end
    end
    % ////////////// MEMORY MANIPULATION (INTERVAL OVERRIDES) /////////////
    methods
        % MEMORY SORTER (INTERVAL)
        function [this] = UpdateMemoryOrderByField(this,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            assert(isfield(this.MEMORY,field),'Field must belong to the memory structure.');
            
            % Handle interval parameters
            if isintval([this.MEMORY.(field)])
                % Reorder the memory structure based on fieldname
                [~,ind] = sort(mid([this.MEMORY.(field)]),2,'descend');     % Ordered indices of the thisect IDs
            else
                [~,ind] = sort([this.MEMORY.(field)],2,'descend');          % Ordered indices of the object IDs
            end
            % Sort the memory structure
            this.MEMORY = this.MEMORY(ind);
        end
        % GET OBJECT PRIORITY (based on distance)
        function [priority_j] = GetObjectPriority(this,objectID)
            % This function calculates a value used to represent the
            % objects priority based on its position. This is called during
            % the object update function.
            logicalIDIndex = [this.MEMORY(:).objectID] == objectID;
            priority_j = 1/mid(this.MEMORY(logicalIDIndex).range(this.MEMORY(logicalIDIndex).sampleNum));
        end
        % GET EMPTY MEMORY STRUCTURE (3D trajectories) [AGENT OVERRIDE]
        function [memStruct]  = CreateMEMORY(this,horizonSteps)
            % This function contains a basic agent-memory structure. This
            % is used to retain information on observed objects and maintain
            % a regular structure.
            
            % Input sanity check
            if nargin < 1
                horizonSteps = 10; % Duration retained in memory
            end
            
            % Handle interval memory types (2D & 3D)
            if this.Is3D
                dim = 3;
            else
                dim = 2;
            end

            % The fields of memory structure define the fields of the
            % simulation 'observation' structure that are retained.
            
            % ///// Create empty memory structure /////
            % House-keeping
            memStruct = struct();
            memStruct.name = ''; 
            memStruct.objectID = uint8(0);
            memStruct.type = OMAS_objectType.misc;
            memStruct.sampleNum = uint8(1);
            % Cartesian measurements
            memStruct.time      = circularBuffer(NaN(1,horizonSteps));
            memStruct.position  = circularBuffer(intval(NaN(dim,horizonSteps)));
            memStruct.velocity  = circularBuffer(intval(NaN(dim,horizonSteps)));
            memStruct.radius    = circularBuffer(intval(NaN(1,horizonSteps)));
            % Spherical measurements
            memStruct.range     = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.heading   = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.elevation = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.width     = circularBuffer(intval(NaN(1,horizonSteps)));
            % Additionals
            memStruct.geometry = struct('vertices',[],'faces',[],'normals',[],'centroid',[]);
            memStruct.priority = [];
            % Trackint (state) containers
%             memStruct.x_k = circularBuffer(intval(NaN(x_length,horizonSteps)));       % Velocity at k
%             memStruct.x_k_plus = circularBuffer(intval(NaN(x_length,horizonSteps)));  % Velocity at k+1 (prediction)
        end
    end
end