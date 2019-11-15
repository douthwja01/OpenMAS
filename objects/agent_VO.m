%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_VO.m) %%%%%%%%%%%%%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_VO < agent
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % AVOIDANCE PARAMETERS
        epsilon = 1E-5;                     % Some reasonably small tolerance
        neighbourDist = 15;                 % The horizon of avoidance consideration
        maxNeighbours = 10;                 % The maximum number of considerable obstacles
        feasabliltyMatrix;                  % The stored achievable velocity matrix     
        pointDensity = 15;                  % Must be odd to allow a zero value
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_VO(varargin)

            % Call the superclass
            this@agent(varargin);  
                        
            % VELOCITY OBSTACLE PARAMETERS
            this.feasabliltyMatrix = this.GetFeasabilityGrid(...
                -ones(3,1)*this.v_max,...
                ones(3,1)*this.v_max,...
                this.pointDensity); 
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
            [this.SENSORS] = this.GetDefaultSensorParameters();              % Default sensing
            %[obj.SENSORS] = obj.GetCustomSensorParameters();              % Experimental sensing
            % /////////////////////////////////////////////////////////////
                             
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup - X = [x;x_dot]' 3D STATE VECTOR
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % BUILD THE STATE VECTOR FOR A 3D SYSTEM WITH CONCATINATED VELOCITIES
            [this] = this.setup_3DVelocities(localXYZVelocity,localXYZrotations);
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
                        
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + this.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            [headingVector] = this.GetTargetHeading();
            desiredVelocity = headingVector*this.v_nominal;

            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet;agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = this.GetAvoidanceCorrection(desiredVelocity,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm      
            
            % /////// COMPUTE STATE CHANGE FROM CONTROL INPUTS ////////////
            [this] = this.Controller(ENV.dt,desiredVelocity);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            this.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = [this.localState(7);this.localState(4:6)];         % Record the control inputs
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    
    % //////////////////// VELOCITY OBSTACLE METHODS //////////////////////
    methods
        % GET THE VO VELOCITY CORRECTION
        function [headingVector,speed] = GetAvoidanceCorrection(this,desiredVelocity,visualiseProblem)
            % This function calculates the collision avoidance heading and
            % velocity correction.
            
            % AGENT KNOWLEDGE
            [p_i,v_i,r_i] = this.GetAgentMeasurements();                    % Its own position, velocity and radius
            
            % Define the obstacle list
            obstacleIDs = [this.MEMORY([this.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];

            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:numel(obstacleIDs)
                % Obstacle data
                p_j = this.GetLastMeasurementByID(obstacleIDs(item),'position');
                
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < this.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < this.neighbourDist;       % [CONFIRMED] 
                if ~neighbourConditionA || ~neighbourConditionB
                    continue
                end
                
                % More obstacle information
                v_j = this.GetLastMeasurementByID(obstacleIDs(item),'velocity');
                r_j = this.GetLastMeasurementByID(obstacleIDs(item),'radius');                                
                % Make none-relative
                p_j = p_j + p_i; 
                v_j = v_j + v_i;                                           % Convert relative parameters to absolute
                tau_b = 0;
                
                % DEFINE THE VELOCITY OBSTACLE PROPERTIES
                VO_i = this.define3DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_b,visualiseProblem);
                % Add a unique identifier
                VO_i.objectID = obstacleIDs(item);             
                VO = vertcat(VO,VO_i);
            end

            % GET THE CAPABLE VELOCITIES (GIVEN AS ABSOLUTE VELOCITIES                
            capableVelocities = this.feasabliltyMatrix;                     % Get all the feasible velocities this timestep 

            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
            escapeVelocities = this.GetEscapeVelocities(capableVelocities,VO); % The viable velocities from the capable 
            
            % APPLY THE MINIMUM DIFFERENCE SEARCH STRATEGY
            [avoidanceVelocity] = this.strategy_minimumDifference(desiredVelocity,escapeVelocities);
                
            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0;0];  % Retain previous heading
            end

            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem && this.objectID == visualiseProblem
                % PLOT THE LEADING TANGENT VECTOR
                OMAS_geometry.drawTriad(p_a,eye(3));
                
                % CURRENT VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),v_a(1),v_a(2),v_a(3),'m');
                q.AutoScaleFactor = 1;
                % DESIRED VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),desiredVelocity(1),desiredVelocity(2),desiredVelocity(3),'g');
                q.AutoScaleFactor = 1;   
                % AVOIDANCE VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),avoidanceVelocity(1),avoidanceVelocity(2),avoidanceVelocity(3),'b');
                q.AutoScaleFactor = 1;
                
                plotCapableVelocities = capableVelocities + p_a;
                plotEscapeVelocities = escapeVelocities + p_a;
                plotAvoidanceVelocities = avoidanceVelocity + p_a;
                scatter3(gca,plotCapableVelocities(1,:),plotCapableVelocities(2,:),plotCapableVelocities(3,:),'r'); 
                scatter3(gca,plotEscapeVelocities(1,:),plotEscapeVelocities(2,:),plotEscapeVelocities(3,:),'g');                 
                scatter3(gca,plotAvoidanceVelocities(1),plotAvoidanceVelocities(2),plotAvoidanceVelocities(3),'b','filled'); 
            end
        end
        
        % GET THE VIABLE ESCAPE VELOCITIES
        function [escapeVelocities] = GetEscapeVelocities(this,velocityMatrix,VOlist)
            % This function takes a matrix of potential velocities and
            % compares them against the agents Velocity Obstacle set. A
            % matrix of escape velocities is then produced.
            % INPUTS:
            % velocityMatrix   - The feasible velocity matrix
            % VOlist           - The velocity obstacle list
            % OUTPUTS:
            % escapeVelocities - The available escape velocity matrix
            
            escapeVelocities = velocityMatrix;  % Prepare output container
            dim = size(escapeVelocities,1);
            
            % FOR EACH VIABLE VELOCITY, CHECK AGAINST CONES
            for i = 1:size(velocityMatrix,2)
                % Nominate a 3D velocity coordinate
                candidatePoint = velocityMatrix(:,i); 
                % Check points against the VO set
                VOnumber = 1; flag = 0;
                while flag ~= 1 && VOnumber <= length(VOlist) 
                    % COLLISION ANNOMOLY
                    if VOlist(VOnumber).openAngle < 0.9*pi
                        % DETERMINE WHETHER POINT IS WITHIN VO
%                        [flag] = obj.isInsideVO(VOlist(VOnumber),candidatePoint);
                        [flag] = this.isInCone(VOlist(VOnumber),candidatePoint);
                    end
                    
                    % Move to the next VO
                    VOnumber = VOnumber + 1;
                end
                
                % If the point lies outside all obstacles %f
                % IF 'flag', point is inside matrix
                if flag 
                   escapeVelocities(:,i) = NaN(dim,1);
                end
            end
            % LOGICALLY REMOVE THE INVALID ELEMENTS
            escapeVelocities = escapeVelocities(:,~isnan(escapeVelocities(1,:)));
        end
        % COMPARE POINT TO CONE GEOMETRY
        function [flag] = isInCone(this,VO,probePoint)      
            % INPUTS:
            % v_b
            % mod_VOaxis
            % unit_VOaxis
            % coneOpenAngle - Cone open angle
            % probePoint    - Test point 
            
            % Compare the angle made between the probeVector and the
            % unitVO axis. Is less than alpha? else not in cone.
                        
            % FROM THE APEX TO THE PROBE POINT
            probeVector = probePoint - VO.apex;                            % Define validation vector
            [unit_probeVector] = probeVector/norm(probeVector);
            
            probeDot = dot(unit_probeVector,VO.axisUnit);
            theta = real(acos(probeDot));                                  % The angle between cone axis and probe vector
                        
            % CHECK POINT RELATION TO CONE GEOMETRY
            flag = 0;
            conditionA = theta - VO.openAngle/2 < this.epsilon;             % With an angle tolerance
            conditionB = probeDot > 0;                                     % In the same direction as the axis vector
            if conditionA && conditionB                                    % Half the cone's open angle
                flag = 1;
            end
        end
        % ASSEMBLE THE 3D VELOCITY OBSTACLE (VO)
        function [VO] = define3DVelocityObstacle(this,p_a,v_a,r_a,p_b,v_b,r_b,tau,plotOn)
            % This function assembles the 3D geometric collision avoidance
            % problem from the relative information observed from the
            % virtual sensor systems.

            % OBSTACLE KNOWLEDGE
            lambda_ab = p_b - p_a;
            r_c = r_a + r_b;        % Configuration radius with a small tolerance offset
            
            % ////////// PARAMETERS DEFINING THE PROBLEM PLANE ////////////
            mod_lambda_ab    = norm(lambda_ab);                            % Scalar relative seperation
            unit_lambda_ab   = lambda_ab/mod_lambda_ab;                    % Unit relative position
            % DEFINE THE ANGLE REFERENCE AXES   
            referenceAxis = [1;0;0]; % To avoid [0;0;0]
            planarNormal = cross(referenceAxis,unit_lambda_ab);             % The reference axis defaults to the agents direction
            if sum(planarNormal) == 0
               planarNormal = cross([0;0;1],unit_lambda_ab);                % Axes then align with the aerospace convention
            end     
            % CALCUATE THE OPEN ANGLE OF THE CONE
            halfAlpha = asin(r_c/mod_lambda_ab);
            halfAlpha = real(halfAlpha);
            
            % CALCULATE THE LEADING TANGENTS ((oldVector,axisVector,theta))
            [leadingTangentVector] = OMAS_geometry.rodriguesRotation(lambda_ab,planarNormal,halfAlpha);

            % CALCULATE THE AXIS PROJECTION
            VOaxis = (dot(leadingTangentVector,lambda_ab)/mod_lambda_ab^2)*lambda_ab; % Define the tangent projection on AB
            axisLength = norm(VOaxis);
            axisUnit = VOaxis/axisLength;
            
            % DEFINE THE LEADING TANGENT VECTOR
            [leadingTangentVector] = OMAS_geometry.rodriguesRotation(VOaxis,planarNormal,halfAlpha);        
            unit_leadingTangent = leadingTangentVector/norm(leadingTangentVector);
            
            % DEFINE THE TRAILING VECTORS (TRAILING)
            [trailingTangentVector] = OMAS_geometry.rodriguesRotation(VOaxis,planarNormal,-halfAlpha);         
            unit_trailingTangent = trailingTangentVector/norm(trailingTangentVector);
            
            % ESTABLISH DIRECTION OF THE AGENT
            leadingTest = unit_leadingTangent'*(v_a - v_b) > unit_trailingTangent'*(v_a - v_b);
%           leadingTest = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))
            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            if leadingTest
                isVaLeading = 1;
            else
                isVaLeading = 0;
            end           
            
            % //////// DEFINE VELOCITY OBSTACLE PARAMETERS ////////////////
            VO = struct('apex',v_b,...
                'axisUnit',axisUnit,...
                'axisLength',mod_lambda_ab,...
                'openAngle',2*halfAlpha,...
                'leadingEdgeUnit',unit_leadingTangent,...
                'trailingEdgeUnit',unit_trailingTangent,...
                'isVaLeading',isVaLeading,...
                'isVaInsideCone',0,...
                'truncationTau',tau,...
                'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
                'truncationCircleRadius',(r_a + r_b)/tau);
            % /////////////////////////////////////////////////////////////
            % IS Va INSIDE THE VO
            [VO.isVaInsideCone] = this.isInCone(VO,v_a); 
            
            % PROBLEM VISUALISATION
            if plotOn && this.objectID == 1               
                % PLOT THE CONE CONSTRUCTION
                coneApex = p_a + VO.apex;
                coneCenter = coneApex + VO.axisUnit*VO.axisLength; 
                % DRAW OBSTACLE AS INPUTTED
                [geometry] = OMAS_graphics.defineSphere(p_b,r_b);
                % REPRESENT GEOMETRY AS A PATCH
                patch(gca,...
                    'Vertices',geometry.vertices,...
                    'Faces',geometry.faces,...
                    'FaceColor','b',...
                    'EdgeColor','k',...
                    'EdgeAlpha',0,...
                    'FaceAlpha',0.2,...
                    'FaceLighting','gouraud',...
                    'LineWidth',0.1);         
                     
                % OBSTACLE VELOCITY (ABS)
                q = quiver3(gca,p_b(1),p_b(2),p_b(3),...
                                v_b(1),v_b(2),v_b(3),'b','filled');          % Plot the obstacles velocity from its relative position
                q.AutoScaleFactor = 1;  
                                     
                % DRAW OBSTACLE AS INPUTTED
                [geometry] = OMAS_graphics.defineSphere(coneCenter,(r_a + r_b));                        
                % REPRESENT GEOMETRY AS A PATCH
                patch(gca,...
                    'Vertices',geometry.vertices,...
                    'Faces',geometry.faces,...
                    'FaceColor','g',...
                    'EdgeColor','k',...
                    'EdgeAlpha',0,...
                    'FaceAlpha',0.2,...
                    'FaceLighting','gouraud',...
                    'LineWidth',0.1);  
                
                % //////////////// VO CONSTRUCTION ////////////////////////
                % THE LEADING TANGENT                
                leadingTangent = VO.leadingEdgeUnit*VO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                            leadingTangent(1),leadingTangent(2),leadingTangent(3),'k');                  
                q.AutoScaleFactor = 1;
                
                % THE TRAILING TANGENT                
                trailingTangent = VO.trailingEdgeUnit*VO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                            trailingTangent(1),trailingTangent(2),trailingTangent(3),'k');
                q.AutoScaleFactor = 1;
                
                % ASSEMBLE VELOCITY OBSTACLE CONES
                try
                    tangentPoint = coneApex + leadingTangentVector;
                    nodes = 11;
                    this.vectorCone(coneApex,coneCenter,tangentPoint,nodes,'g');
                catch
                end
            end
        end
    end   
    % ///////////////// STATIC VELOCITY OBSTACLE METHODS //////////////////
    methods (Static)
      % DECOMPOSE THE OBSTACLE GEOMETRIES INTO A VO SET
        function [VO] = define3DComplexVelocityObstacle(p_a,v_a,r_a,p_b,v_b,geometry,tau)
            % This function computes the VO set from a list of complex
            % obstacles defined by their .GEOMETRY parameter.
            
            % NOTES:
            % - The .geometry parameter is the geometry of the object pre-
            %   positioned relative to the agent, rotated and scaled.
            % - VO structure:
            %   VO = struct('apex',v_b,...
            %     'axisUnit',axisUnit,...
            %     'axisLength',mod_lambda_ab,...
            %     'openAngle',2*halfAlpha,...
            %     'leadingEdgeUnit',unit_leadingTangent,...
            %     'trailingEdgeUnit',unit_trailingTangent,...
            %     'isVaLeading',isVaLeading,...
            %     'isVaInsideCone',0,...
            %     'truncationTau',tau,...
            %     'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
            %     'truncationCircleRadius',(r_a + r_b)/tau);
            
            norm_lateralPositionProjection = norm([p_b(1:2);0]);
            unit_lateralPositionProjection = [p_b(1:2);0]/norm_lateralPositionProjection;   % XY PLANAR PROJECTION OF RELATIVE POSITION
            
            a_min = 0; a_max = 0;
            % Get the points with the largest perpendicular projection
            for v = 1:size(geometry.vertices,1)
                % GET THE XY VERTEX PROJECTION
                lateralVertexProjection = [geometry.vertices(v,1:2)';0];          % Projection of the vertex in the XY plane
                unit_lateralVertexProjection = lateralVertexProjection/norm(lateralVertexProjection);
                % GET THE ANGULAR ARGUEMENT WITH THE POSITION VECTOR
                crossProduct = cross(unit_lateralVertexProjection,unit_lateralPositionProjection); % Planar normal
                dotProduct   = dot(unit_lateralVertexProjection,unit_lateralPositionProjection);
                % THE SIGNED ANGLE OF THE VERTEX RELATIVE TO THE
                % CENTROID VECTOR
                signedAngularProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
                % RETAIN THE LARGEST HEADING CHANGE AS THE VO BOUNDARY
                if signedAngularProjection > a_max
                    a_max = signedAngularProjection;
                    unit_trailingTangent = unit_lateralVertexProjection; % Assuming that clockwise +ve
                    %                        maximalVertexProjection = [obstacleGeometry.vertices(v,1:2)';0];
                elseif signedAngularProjection < a_min
                    a_min = signedAngularProjection;
                    unit_leadingTangent = unit_lateralVertexProjection;  % Assuming that anti-clockwise -ve
                    %                        minimalVertexProjection = [obstacleGeometry.vertices(v,1:2)';0];
                end
            end
            
            % CHECK IF Va belongs to the projection cone
            crossProduct = cross(v_a,unit_lateralPositionProjection);      % Planar normal
            dotProduct = dot(v_a/norm(v_a),unit_lateralPositionProjection);% Angle between two
            signedVelocityProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
            % DETERMINE IF THE CURRENT VELOCITY BELONGS TO THE CURRENT VO
            isVaInsideCone = 0;
            if signedVelocityProjection < a_max && signedVelocityProjection > a_min
                isVaInsideCone = 1;
            end
            
            % THE MAXIMAL VERTEX PROJECTIONS
            equivalentOpenAngle = abs(a_min) + abs(a_max);
            isVaLeading = 1;                                               % Va always leading (obstacle is static)
            
            effectiveRadii = norm_lateralPositionProjection*sin(equivalentOpenAngle);       % <<< TO BE DEFINED         %(r_a + r_b)
            % DEFINE THE VO STRUCTURE
            VO = struct('apex',v_b,...
                'axisUnit',unit_lateralPositionProjection,...
                'axisLength',norm_lateralPositionProjection,...
                'openAngle',equivalentOpenAngle,...
                'leadingEdgeUnit',unit_leadingTangent,...
                'trailingEdgeUnit',unit_trailingTangent,...
                'isVaLeading',isVaLeading,...
                'isVaInsideCone',isVaInsideCone,...
                'truncationTau',tau,...
                'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
                'truncationCircleRadius',(effectiveRadii)/tau);            % Double the radii of A clear
        end
        % SMALLEST DEVIATION FROM DESIRED 
        function [optimalVelocity] = strategy_minimumDifference(desiredVelocity,escapeVelocities)
            % This function searches a matrix of viable velocities by 
            % applying the minimum difference from the desired trajectory 
            % approach.  
            
            % ESTABLISH THE INPUT DIMENSIONS
            inputDim = numel(desiredVelocity);        
            
            % CONSTRUCT THE SEARCH SPACE
            searchMatrix = vertcat(escapeVelocities,zeros(1,size(escapeVelocities,2)));
            
            for i = 1:size(searchMatrix,2)
                searchMatrix((inputDim+1),i) = norm(desiredVelocity - searchMatrix(1:inputDim,i)); % Add a dimension with the absolute vector value
            end
            % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
            [~,minIndex] = min(searchMatrix((inputDim+1),:),[],2);                    % Return the index of the smallest vector
            optimalVelocity = searchMatrix(1:inputDim,minIndex);
            
            % ERROR CATCHING
            if isempty(optimalVelocity)
                warning('No viable velocities found in search matrix.');
                optimalVelocity = zeros(inputDim,1);
            end            
        end
        % BUILD A UNIT, CUBOID N-DIMENSIONAL POINT CLOUD
        function [cubePoints] = GetFeasabilityGrid(minVector,maxVector,pointDensity)
            % This function generates an n-dimensional cube of points defined
            % by the minimum and maximum vector arguements. 
            % INPUTS:
            % minVector    - A vector of minimum bounds
            % maxVector    - A vector of maximum bounds
            % pointDensity - The number of descrete steps 
            % OUTPUTS:
            % cubePoints   - A vector of the 3D points within the cloud
            
            % CONSTANTS
            dimensionality = numel(minVector); 
            numPoints = pointDensity^dimensionality;
            % CONTAINERS
            dimMultipliers = zeros(dimensionality,1);
            dimensionalDistribution = zeros(dimensionality,pointDensity);
            for dim = 1:dimensionality
                dimensionalDistribution(dim,:) = linspace(minVector(dim),maxVector(dim),pointDensity);  % Matrix with the dimensional distribution
                dimMultipliers(dim) = numPoints/pointDensity^dim;                                       % Dimension tessilation indicies
            end     
            dimMultipliers = fliplr(dimMultipliers);
            
            % GENERATE THE N-DIMENSIONAL POINT LIST
            cubePoints = zeros(dimensionality,size(dimensionalDistribution,2)^dimensionality);
            for dim = 1:dimensionality
                insertVector = [];
                for i = 1:pointDensity
                   distributionIter = repmat(dimensionalDistribution(dim,i),1,dimMultipliers(dim));
                   insertVector = horzcat(insertVector,distributionIter);
                end
                no_copies = numPoints/size(insertVector,2);
                cubePoints(dim,:) = repmat(insertVector,1,no_copies);
            end 
        end
    end
    % ///////////////////////////// UTILITIES /////////////////////////////
    methods (Static)
        % DRAW CONE
        function [Cone]  = vectorCone(pointA,pointB,radialPoint,nodes,coneColour)
            % This function is designed to construct a cone mesh between two 
            % 3D points Using a radial Point to specify the radial properties.
            % INPUTS:
            % pointA - First cone axis point (origin).
            % pointB - Second cone axis point.
            % radialPoint - Third point in space that defines the cones maximum radius
            % nodes  - The number of mesh nodes
            % OUTPUTS
            % Cone   - Cone mesh structure
            
            % INPUT HANDLING
            if ~exist('nodes','var')
                nodes = 10;
            end
            
            % MANUAL INPUTS
            coneEdgeColour = 'r';
            coneAlpha = 0.1;
            coneClosed = 0;
            coneLines = 0;
            
            % DEFINE VECTOR PROBLEM PARAMETERS
            axisVector = pointB - pointA;
            mod_axisVector = sqrt(sum((axisVector).^2)); % Axis vector properties
            tangent = radialPoint - pointA;
            mod_tangent = sqrt(sum((tangent).^2));       % Tangental vector properties
            
            % DEFINE THE TANGENT-AXIS PROJECTION
            trueAB = (dot(tangent,axisVector)/mod_axisVector^2)*axisVector;
            mod_trueAB = sqrt(sum((trueAB).^2));
            
            % GET THE RADIUS MODULUS
            mod_radius = sqrt(mod_tangent^2-mod_trueAB^2);
                                
            % Creating 2 circles in the YZ plane
            t=linspace(0,2*pi,nodes)';      % Create axis point set
            xa2 = zeros(length(t),1); 
            xa3 = zeros(size(xa2)); 
            xb2 = mod_radius*cos(t);
            xb3 = mod_radius*sin(t);        % Scale the second axis in the 
            
            % Creating the points in the X-Direction
            x1=[0 mod_trueAB];
            
            % Creating (Extruding) the cylinder points in the X-Directions
            xx1 = repmat(x1,length(xa2),1);
            xx2 = [xa2 xb2];  
            xx3 = [xa3 xb3]; % Concatinate second circle set
            
            % Drawing two filled cirlces to close the cylinder
            if coneClosed == 1
                EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
                EndPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'r');
            end
            
            % GENERATE THE CONE 
            % Plot the cone from the origin along x-axis and scale to size 
            Cone = mesh(gca,real(xx1),real(xx2),real(xx3));
            
            % Get the planar rotation angle
            unit_Vx=[1 0 0];
            angle_X1X2 = acos(dot(unit_Vx,axisVector)/(norm(unit_Vx)*mod_axisVector))*180/pi;
            
            % Get rotation axis
            axis_rot = cross([1 0 0],axisVector); 
            
            % Rotating the plotted Cone and the end plate circles to the required
            % angles
            if angle_X1X2~=0 % Rotation is not needed if required direction is along X
                rotate(Cone,axis_rot,angle_X1X2,[0 0 0])
                if coneClosed==1
                    rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
                    rotate(EndPlate2,axis_rot,angle_X1X2,[0 0 0])
                end
            end
            
            % Till now Cone has only been aligned with the required direction, but
            % position starts from the origin. so it will now be shifted to the right
            % position
            if coneClosed == 1
                set(EndPlate1,'XData',get(EndPlate1,'XData') + pointA(1))
                set(EndPlate1,'YData',get(EndPlate1,'YData') + pointA(2))
                set(EndPlate1,'ZData',get(EndPlate1,'ZData') + pointA(3))
                
                set(EndPlate2,'XData',get(EndPlate2,'XData') + pointA(1))
                set(EndPlate2,'YData',get(EndPlate2,'YData') + pointA(2))
                set(EndPlate2,'ZData',get(EndPlate2,'ZData') + pointA(3))
            end
            set(Cone,'XData',get(Cone,'XData') + pointA(1))
            set(Cone,'YData',get(Cone,'YData') + pointA(2))
            set(Cone,'ZData',get(Cone,'ZData') + pointA(3))
            
            % SET THE COLOUR OF THE CONE AND END PLATES
            set(Cone,'AmbientStrength',1,...
                     'FaceColor',coneColour,...
                     'FaceLighting','gouraud',...
                     'FaceAlpha',coneAlpha,...
                     'EdgeColor',coneEdgeColour);        % Cone verticies
            if coneClosed==1
                set([EndPlate1 EndPlate2],...
                    'AmbientStrength',1,...
                    'FaceColor',coneColour,...
                    'FaceLighting','gouraud',...
                    'FaceAlpha',coneAlpha);         % End-plate 
            end
            
            % If lines are not needed making it disapear
            if coneLines == 0
                set(Cone,'EdgeAlpha',0)
            end
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]