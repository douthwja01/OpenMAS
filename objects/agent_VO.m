%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_VO.m) %%%%%%%%%%%%%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_VO < agent
    
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES OF THE VELOCITY OBSTACLE METHOD
        feasabliltyMatrix;                  % The stored achievable velocity matrix
        pointDensity;                       % Must be odd to allow a zero value
        obstacleSafetyFactor;               % Modify the apparent size of the obstacle
        % gif toggle
        giffOn = 0;
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_VO(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'
            % CONTROL PARAMETERS
            obj.priorError = zeros(4,1);
            % PARAMETERISE KINEMATIC CONSTRAINTS
            obj.linearVelocityLimits = [15;15;15];                         % Limits on the agents velocity 
            obj.linearAccelerationLimits = [10;10;10];                     % Limits on the agents acceleration
            obj.angularVelocityLimits = [inf;inf;inf];
            obj.angularAccelerationLimits = [inf;inf;inf]; % rad/s^2
            
            % VELOCITY OBSTACLE PARAMETERS  
            obj.obstacleSafetyFactor = 1.2;                                % Modify the apparent size of the obstacle       
            obj.pointDensity = 11;                                          % Must be odd to allow a zero value
            obj.feasabliltyMatrix = obj.getFeasabilityGrid(-obj.linearAccelerationLimits,...
                                                            obj.linearAccelerationLimits,...
                                                            obj.pointDensity);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        
        %% AGENT MAIN CYCLES
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            observationSet = varargin{1};                                  % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(observationSet);
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredSpeed = 2;
            if ~isempty(obj.targetWaypoint)
                % DEFAULT THE HEADING VECTOR TOWARDS THE WAYPOINT
                waypointPosition = obj.targetWaypoint.state(1:3);
                desiredVelocity  = (waypointPosition/norm(waypointPosition))*desiredSpeed;
            else
                desiredVelocity = [1;0;0]*desiredSpeed;
            end
                        
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % /////////////// PARSE CONTROLLER INPUTS /////////////////////
            targetSpeed = norm(desiredVelocity);
            targetHeading = desiredVelocity/targetSpeed;                   % Get the control inputs
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            useController = 0;
            if useController
                [d_heading,d_speed,obj] = obj.trajectoryController(targetHeading,targetSpeed);  
            else
                [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);
            end
            newSpeed = (norm(obj.localState(7:9)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = d_heading;                                        % Heading change relative to straight ahead
            
            % //////////////// APPLY KINEMATIC LIMITS /////////////////////
            [newHeading,newSpeed] = obj.kinematicContraints(dt,newHeading,newSpeed);
            newVelocity = [1;0;0]*newSpeed;                                % Convert to local velocity vector
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            [newState] = obj.stateDynamics_velocities(dt,newVelocity,newHeading);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                newState = obj.freezeAgent();
            end
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:length(obj.priorError),TIME.currentStep) = [newSpeed;newState(4:6)];         % Record the control inputs 
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
        end
        % GET THE VO VELOCITY CORRECTION
        function [avoidance_V] = getAvoidanceCorrection(obj,dt,desired_V,knownObstacles,visualiseProblem)
            % This function calculates the collision avoidance heading and
            % velocity correction.
            
            % AGENT KNOWLEDGE
            p_a = obj.localState(1:3,1);
            v_a = obj.localState(7:9,1);
            r_a = obj.VIRTUAL.size;
            capable_V = v_a + obj.feasabliltyMatrix*dt;                    % Get all the feasible velocities this timestep 
            
            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:length(knownObstacles)
                % OBSTACLE KNOWLEDGE
                p_b = knownObstacles(item).state(1:3) + p_a; 
                v_b = knownObstacles(item).state(4:6) + v_a;               % Convert relative parameters to absolute
                r_b = knownObstacles(item).size;
                % DEFINE THE VELOCITY OBSTACLE PROPERTIES
                [VOorigin,VOaxis,VOangle] = obj.defineVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem);
                VO = vertcat(VO,struct('origin',VOorigin,'axis',VOaxis,'angle',VOangle));
            end
            
            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
            escape_V = obj.getEscapeVelocities(capable_V,VO);              % The viable velocities from the capable 
            if isempty(escape_V)
                warning('No viable escape velocities');
                avoidance_V = [1;0;0]*norm(desired_V);
                return
            end
            
            % SEARCH THE VIABLE ESCAPE VELOCITIES FOR THE VELOCITY WITH THE 
            % SMALLEST DEVIATION FROM THE DESIRED VELOCITY          
            searchMatrix = vertcat(escape_V,zeros(1,size(escape_V,2)));
            for i = 1:size(searchMatrix,2)
                searchMatrix(4,i) = norm(desired_V - searchMatrix(1:3,i)); % Add a dimension with the absolute vector value
            end
            % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
            [~,minIndex] = min(searchMatrix(4,:),[],2);                    % Return the index of the smallest vector
            avoidance_V = searchMatrix(1:3,minIndex);
            
            % ERROR CATCHING
            if isempty(avoidance_V)
                warning('No avoidance velocities available for agent %s',obj.name);
                avoidance_V = [1;0;0]*norm(desired_V);
            end
               
            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem
                hold on;
                scatter3(capable_V(1,:),capable_V(2,:),capable_V(3,:),'r'); 
                scatter3(escape_V(1,:),escape_V(2,:),escape_V(3,:),'g'); 
                q = quiver3(0,0,0,current_V(1),current_V(2),current_V(3),'b');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,desired_V(1),desired_V(2),desired_V(3),'g');
                q.AutoScaleFactor = 1;               
                q = quiver3(0,0,0,avoidance_V(1),avoidance_V(2),avoidance_V(3),'r');
                q.AutoScaleFactor = 1;
                scatter3(avoidance_V(1,:),avoidance_V(2,:),avoidance_V(3,:),'b','filled'); 
            end
        end
        
        %% AVOIDANCE VELOCITY PARSING
        % GET THE VIABLE ESCAPE VELOCITIES
        function [escapeVelocities] = getEscapeVelocities(obj,velocityMatrix,VOlist)
            % This function takes a matrix of potential velocities and
            % compares them against the agents Velocity Obstacle set. A
            % matrix of escape velocities is then produced.
            % INPUTS:
            % velocityMatrix   - The feasible velocity matrix
            % VOlist           - The velocity obstacle list
            % OUTPUTS:
            % escapeVelocities - The available escape velocity matrix
            
            escapeVelocities = [];   % Prepare output container
            
            % FOR EACH VIABLE VELOCITY, CHECK AGAINST CONES
            for i = 1:size(velocityMatrix,2)
                % Nominate a 3D velocity coordinate
                candidatePoint = velocityMatrix(:,i); 
                % Check points against the VO set
                VOnumber = 1; flag = 0;
                while flag ~= 1 && VOnumber <= length(VOlist)   
                    % Calculate axis properties
                    mod_coneAxis = norm(VOlist(VOnumber).axis);
                    unit_coneAxis = VOlist(VOnumber).axis/mod_coneAxis;
                    
                    % DETERMINE WHETHER POINT IS WITHIN VO
                    [flag] = obj.isInCone(VOlist(VOnumber).origin,unit_coneAxis,...
                                          mod_coneAxis,VOlist(VOnumber).angle,...
                                          candidatePoint);
                    % Move to the next VO
                    VOnumber = VOnumber + 1;
                end
                % If the point lies outside all obstacles
                if ~flag 
                    escapeVelocities = horzcat(escapeVelocities,candidatePoint);
                end
            end
        end
        
        %% VELOCITY OBSTACLE DESCRIPTION METHODS
        % DRAW THE VELOCITY OBSTACLE PROBLEMS
        function [VOorigin,VOaxis,VOgamma] = defineVelocityObstacle(obj,p_a,v_a,r_a,p_b,v_b,r_b,plotOn)
            % This function assembles the geometric collision avoidance
            % problem from the relative information observed from the
            % virtual sensor systems.

            % OBSTACLE KNOWLEDGE
            lambda_ab = p_b - p_a;
            c_b       = v_b - v_a;
            r_c = r_a + obj.obstacleSafetyFactor*r_b;
            
            % ////////// PARAMETERS DEFINING THE PROBLEM PLANE ////////////
            mod_lambda_ab    = norm(lambda_ab);                            % Scalar relative seperation
            unit_lambda_ab   = obj.unit(lambda_ab);                        % Unit relative position
            unit_referenceAxis  = [0;0;1];                                 % Unit reference axis
            % DEFINE THE ANGLE REFERENCE AXES                            
            referenceVector = [1;0;0];                                     % The reference axis defaults to the agents direction
            problemAxis = cross(referenceVector,unit_lambda_ab);
            if problemAxis == 0
               referenceVector = unit_referenceAxis;                       % Else form a comparison to the unit y axis
               problemAxis = cross(referenceVector,unit_lambda_ab);           % Axes then align with the aerospace convention
            end     
                 
            % The projection is defined as the obstacle (in the
            % configuration space) augmented by its relative velocity.
            radiusVO = -r_c*unit_lambda_ab;                   % Unit radius of the obstacle (direction of lambdaAB)
            VOgamma = asin(r_c/mod_lambda_ab);               % Get rotation angle for tangential vector          
                                    
            % DEFINE THE TANGENT VECTORS
            angle = (pi/2) - VOgamma;                                                 % Define the leading tangent angle 
            [leadingTangentialRadius] = obj.rotateVectorAboutAxis(radiusVO,problemAxis,angle); % The leading tangential radial vector         
            leadingTangentVector = lambda_ab + leadingTangentialRadius;      % The vector to the tangent point
            VOaxis = (dot(leadingTangentVector,lambda_ab)/mod_lambda_ab^2)*lambda_ab; % Define the tangent projection on AB
            
            % VELOCITY OBSTACLE VECTOR DESCRIPTION
            VOorigin = v_b;
            VOlambda = VOorigin + lambda_ab;
            VOtangentPoint = VOorigin + lambda_ab + leadingTangentialRadius;
            
            % PROBLEM VISUALISATION
            if plotOn
                % DRAW AGENT FORWARD VELOCITY
                coneCenter = VOorigin + VOaxis;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(0,0,0,1,0,0,'r');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,-1,0,'g');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,0,-1,'b');
                q.AutoScaleFactor = 1;
                
                % DRAW OBSTACLE REPRESENTATION
                [Xb,Yb,Zb] = obj.getSphere(lambda_ab,r_b);                            
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                         'FaceColor','g',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0); 
                % PLOT THE RELATIVE VELOCITY OF THE OBSTACLE
                q = quiver3(lambda_ab(1),lambda_ab(2),lambda_ab(3),...
                            v_b(1), v_b(2), v_b(3),'b','filled');          % Plot the obstacles velocity from its relative position
                q.AutoScaleFactor = 1;                     
                     
                % \\\\\\  DRAW THE VELOCITY OBSTACLE PROJECTION \\\\\\\\\\\
                % PLOT THE VO ORIGIN
                q = quiver3(p_a(1),p_a(2),p_a(3),...
                            VOorigin(1),VOorigin(2),VOorigin(3),'g');
                q.AutoScaleFactor = 1;
                % DRAW THE VELOCITY OBSTACLE PROJECTION        
                [X_vo,Y_vo,Z_vo] = obj.getSphere(VOlambda,r_c);  
                sphB = mesh(X_vo,Y_vo,Z_vo); 
                set(sphB,'facealpha',0.2,...
                         'FaceColor','g',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0);
                % PLOT THE LEADING TANGENTIAL RADIUS OF THE VO
                q = quiver3(VOlambda(1),VOlambda(2),VOlambda(3),...
                            leadingTangentialRadius(1),leadingTangentialRadius(2),leadingTangentialRadius(3),'k');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(VOorigin(1),VOorigin(2),VOorigin(3),...
                            leadingTangentVector(1),leadingTangentVector(2),leadingTangentVector(3),'k');
                q.AutoScaleFactor = 1;
                % ASSEMBLE VELOCITY OBSTACLE CONES
                nodes = 11;
                obj.vectorCone(VOorigin,coneCenter,VOtangentPoint,nodes,'g');
            end
        end               
    end
    % STATIC VELOCITY OBSTACLE METHODS
    methods (Static)
        % BUILD A UNIT, CUBOID N-DIMENSIONAL POINT CLOUD
        function [cubePoints] = getFeasabilityGrid(minVector,maxVector,pointDensity)
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
        % COMPARE POINT TO CONE GEOMETRY
        function [flag] = isInCone(origin,unit_coneAxis,mod_coneAxis,coneAngle,point)      
            % INPUTS:
            % VOorigin
            % mod_VOaxis
            % unit_VOaxis
            % alpha      - Cone angle
            % probePoint - Test point 
            
            % Compare the angle made between the probeVector and the
            % unitVO axis. Is less than alpha? else not in cone.
            
            if numel(origin) ~= numel(point)
                error('dimensionality problem');
            end                
            
            probeVector = point - origin;           % Define validation vector
            [unit_probeVector] = probeVector/sqrt(sum(probeVector.^2));
            theta = acos(dot(unit_probeVector,unit_coneAxis));      % The angle between cone axis and probe vector
            
            % CHECK POINT RELATION TO CONE GEOMETRY
            flag = 0;
            if theta <= coneAngle
                probeVOProjection = dot(probeVector,unit_coneAxis); %Calculate the projection of the probeVector onto VOaxis
                if probeVOProjection <= mod_coneAxis
                    % Is the projection of the point less than the cone length
                    % YES, IN CONE
                    flag = 1;
                end
            end
        end
        % DRAW CONE
        function [Cone] = vectorCone(pointA,pointB,radialPoint,nodes,coneColour)
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
            Cone = mesh(real(xx1),real(xx2),real(xx3));
            
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
        % DRAW SPHERE
        function [X,Y,Z] = getSphere(position,radius)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius
            [X,Y,Z] = sphere(40);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
        end
        % CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
        function [outputMatrix] = skew(inputVector)
            % This function generates a skew-symmetric for the computation
            % of the vector cross-product.
            % INPUTS:
            % inputVector - The original 3D vector
            % OUTPUT:
            % outputMatrix - The equivalent skew-symmetric matrix
            
            if length(inputVector) ~= 3
                warning('The input vector must be three dimensional.');
                return
            end
            % Apply element mapping
            outputMatrix = zeros(3,3);
            outputMatrix(1,2) = -inputVector(3);
            outputMatrix(1,3) =  inputVector(2);
            outputMatrix(2,1) =  inputVector(3);
            outputMatrix(2,3) = -inputVector(1);
            outputMatrix(3,1) = -inputVector(2);
            outputMatrix(3,2) =  inputVector(1); % Arrange the components
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]