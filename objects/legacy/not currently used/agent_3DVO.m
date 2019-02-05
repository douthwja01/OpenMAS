%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_3DVO.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_3DVO < agent
    
    %% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % AVOIDANCE TERMS
        obstacleSafetyFactor;                                              % Modify the apparent size of the obstacle
        % CONTROLLER TERMS
        priorError = zeros(4,1);
        % gif toggle
        giffOn = 0;
    end
    %  CLASS METHODS
    methods
        % CONSTRUCTOR
        function obj = agent_3DVO(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);     
            % PARAMETERISE KINEMATIC CONSTRAINTS
            obj.linearVelocityLimits = [30;30;1];                          % Override the 'agent' linear velocity constraints 
            obj.linearAccelerationLimits = [1;1;1];                        % Override the 'agent' linear acceleration constraints
            obj.angularVelocityLimits = [pi/3;pi/3;pi/3];                  % Override the 'agent' angular velocity constraints
            obj.angularAccelerationLimits = [0.2;0.2;0.2];                 % Override the 'agent' angular acceleration constraints
            % AGENT SPECIFIC PARAMETERS
            obj.obstacleSafetyFactor = 1.0;                                % Modify the apparent size of the obstacle
            obj.sensorRange = 50;                                         % Virtual sensor range
            % VIRTUAL DEFINITION
            obj.VIRTUAL.detectionRange = obj.sensorRange;                  % Assign the range attribute to the SIM VIRTUAL attribute
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);                     % Parse inputs for specific parameter overrides
        end
        % AGENT MAIN CYCLE
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % \\\\\\\\\\\\\\\\\\ ASSIGN DEFAULT BEHAVIOUR \\\\\\\\\\\\\\\\\
            desiredSpeed = 12;                                             % 20m/s
            
            % CHECK FOR NEW OBSTACLES & UPDATE WAYPOINT KNOWLEDGE
            observationSet = varargin{1};                                  % The detected objects
            [obj,obstacleSet,~] = obj.getAgentUpdate(observationSet);      % Update agents 'memory' structure
            
            % GET WAYPOINT (TARGET) HEADING VECTOR 
            if ~isempty(obj.targetWaypoint)  
               % GET THE WAYPOINT HEADING VECTOR
               desiredHeading = obj.targetWaypoint.state(1:3)/norm(obj.targetWaypoint.state(1:3));
            else
               desiredHeading = [1;0;0];                                   % Directly forward
            end
            
            % \\\\\\\\\\\\\\\\\\\\\ FIGURE CONTAINER \\\\\\\\\\\\\\\\\\\\\\
            visualiseAgent = 1; visualiseProblem = 0;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis vis3d;
                view([-78 50]);
                xlabel('x_{m}');
                ylabel('y_{m}');
                zlabel('z_{m}');
            end
            % \\\\\\\\\\\\\\\\\\\\\ OBSTACLE AVOIDANCE \\\\\\\\\\\\\\\\\\\\
            algorithm_start = tic; algorithm_indicator = 0;                % Avoidance routine not yet ran
            if ~isempty(obstacleSet)                                 
                % IF THE OBSTACLE IS A PROBLEM
                algorithm_indicator = 1;
                for item = 1:length(obstacleSet)
                    % GET THE OBSTACLES INFORMATION
                    sizeEstimateB     = obstacleSet(item).size;            % Get object size from memory
                    positionEstimateB = obstacleSet(item).state(1:3);      % Get current position from memory
                    velocityEstimateB = obstacleSet(item).state(4:6);      % Get current velocity from memory
                    % GET COLLISION LIKELYHOOD
                    [tau] = obj.validateCollision(positionEstimateB,velocityEstimateB);
                    if any(isnan(velocityEstimateB)) || tau < 0
                       continue                                       % Avoidance requires a velocity estimate
                    end
                    
                    
                end
                
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the avoidance routine


            % \\\\\\\\\\\\\\ CALCULATE CONTROL INPUTS \\\\\\\\\\\\\\\\\\\\\
            % This controller is used to drive an error in forward speed
            % and heading to zero.
            
            % DEFINE THE CONTROL ERROR [e_x;e_phi;e_theta;e_psi]
            speedError = desiredSpeed - obj.localState(7);                 % Get the linear (x) speed error
            [thetaError,lambdaError] = obj.getControlInputs([1;0;0],desiredHeading); % Get the angular heading error
            controlError = [speedError;0;thetaError;-lambdaError];         % -ve in the NED control frame of reference
            % BUILD THE CONTROL INPUTS
            Kp_lin = 3;   Kd_lin = 0.1;
            Kp_ang = 0.3; Kd_ang = 10;  
            control_inputs = diag([Kp_lin;0;Kp_ang;Kp_ang])*controlError + ...
                             diag([Kd_lin;0;Kd_ang;Kd_ang])*(controlError-obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = controlError;

            % \\\\\\\\\\\\\\ APPLY KINEMATIC LIMITS \\\\\\\\\\\\\\\\\\\\\\\
            [control_inputs(1)] = obj.boundValue(control_inputs(1),-obj.linearAccelerationLimits(1),obj.linearAccelerationLimits(1));
            [control_inputs(2)] = obj.boundValue(control_inputs(2),-obj.angularAccelerationLimits(1),obj.angularAccelerationLimits(1));
            [control_inputs(3)] = obj.boundValue(control_inputs(3),-obj.angularAccelerationLimits(2),obj.angularAccelerationLimits(2));
            [control_inputs(4)] = obj.boundValue(control_inputs(4),-obj.angularAccelerationLimits(3),obj.angularAccelerationLimits(3));

            % 7. PLOT THE AVOIDANCE SITUATION
            if visualiseProblem && obj.objectID == visualiseAgent
                % DRAW AGENT SPHERE
                [Xb,Yb,Zb] = obj.getSphere([0;0;0],obj.VIRTUAL.size);
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                    'FaceColor','b',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % DRAW LOCAL TRIAD (FOR REFERENCE)
                q = quiver3(0,0,0,1,0,0,'r');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,-1,0,'g');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,0,-1,'b');
                q.AutoScaleFactor = 1;
                
                %  PLOT THE LOCAL AGENT VELOCITY
                q = quiver3(0,0,0,obj.localState(7),obj.localState(8),obj.localState(9),'b');   % Plot the current velocity
                q.AutoScaleFactor = 1;
                
                %  PLOT THE DESIRED AGENT VELOCITY
                desiredVelocity = desiredSpeed*desiredHeading;
                q = quiver3(0,0,0,desiredVelocity(1),desiredVelocity(2),desiredVelocity(3),'r');
                q.AutoScaleFactor = 1;
                view([-45 45 ]);
                drawnow;
                % SEND THE FIGURE TO THE GIF GENERATOR
                [obj] = obj.getAnimationFrame(overHandle,TIME,'3DVO_agent');
                close(overHandle);
            end
            
            % \\\\\\\\\\\\\\\ GET THE NEW STATE VECTOR \\\\\\\\\\\\\\\\\\\\
            linear_inputs  = [1;0;0]*control_inputs(1);
            angular_inputs = control_inputs(2:4);
            % \\\\\\\\\\\\\\\ UPDATE THE DYNAMICS     
            newState = obj.stateDynamics(dt,linear_inputs,angular_inputs);
            
            % \\\\\\\\\\\\\\\ UPDATE THE CLASS GLOBAL PROPERTIES \\\\\\\\\\
            obj = obj.updateGlobalProperties(dt,newState);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:4,TIME.currentStep) = control_inputs;               % Record the control inputs
            obj.DATA.error(1:4,TIME.currentStep) = controlError;
        end
    end
    
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods
        % DRAW THE VELOCITY OBSTACLE PROBLEMS
        function [AVO,DVO,alphaVO]  = defineVelocityObstacle(obj,obstacle,plotOn)
            % This function assembles the geometric collision avoidance
            % problem from the relative information observed from the
            % virtual sensor systems.
            % OUTPUTS:
            % AVO - The VO origin coordinates
            % DVO - The projection of the tangent on Doi (vector to the center of the plane containing all tangents)
            % alphaVO - The cone projection angle (rad)
            
            % DECLARE STATIC OBSTACLE KNOWLEDGE
            % Obstacle information
            Doi  = obstacle.state(1:3);                                    % Relative position (center to center)
            Vr   = obstacle.state(4:6);                                    % Relative velocity
            mDoi = sqrt(sum(Doi.^2));                                      % Scalar relative seperation
            UDoi = obj.unit(Doi);                                          % Relative Heading
            
            % DEFINE THE ANGLE REFERENCE AXES
            referenceHeading = [1;0;0];                                    % The reference axis defaults to the agents direction
            planarAxis = cross(referenceHeading,UDoi);
            if planarAxis == 0
                planarAxis = cross([0;0;1],UDoi);                          % Get the vector defining the plane the problem will be solved in
            end
            
            % DETERMINE THE VELOCITY OBSTACLE PROJECTION
            origin = [0;0;0];
            Dvo = origin + Doi + Vr;                                       % The vector to the center point of the translated sphere
            
            % PROTECTED ZONE RADIUS
            rpz = obj.VIRTUAL.size + obstacle.size*obj.obstacleSafetyFactor;

            Rvo = -rpz*UDoi;                                               % Unit radius of the obstacle (direction of doi)
            alphaVO = asin(rpz/mDoi);                                      % Get rotation angle for tangential vector
            
            % DEFINE THE TANGENT VECTORS
            angle = (pi/2) - alphaVO;                                      % Define the leading tangent angle
            [tangRadiusL] = obj.rotateVectorAboutAxis(Rvo,planarAxis,angle);
            LT = Doi + tangRadiusL;                                        % The leading tangent vector
            LTP = Dvo + tangRadiusL;                                       % The leading tangent point
            
            %% DEFINE FINAL VELOCITY OBSTACLE ATTRIBUTES
            AVO = Vr;                                                      % Define the VO origin point
            DVO = (dot(LT,Doi)/mDoi^2)*Doi;                                % Define the tangent projection on AB
            
            % PROBLEM VISUALISATION
            if plotOn               
                % DRAW OBSTACLE REPRESENTATION
                [Xb,Yb,Zb] = obj.getSphere(Doi,obstacle.size);
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                    'FaceColor','g',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % DRAW THE PROJECTED OBSTACLE 
                [X_vo,Y_vo,Z_vo] = obj.getSphere(Dvo,rpz);
                sphB = mesh(X_vo,Y_vo,Z_vo);
                set(sphB,'facealpha',0.2,...
                    'FaceColor','r',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % PLOT THE RELATIVE VELOCITY OF THE OBSTACLE
                q = quiver3(Doi(1),Doi(2),Doi(3),...
                             Vr(1), Vr(2), Vr(3),...   % Plot the obstacles velocity from its relative position
                               'b','filled');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENTIAL RADIUS OF THE VO
                q = quiver3(Dvo(1),Dvo(2),Dvo(3),...
                            tangRadiusL(1),tangRadiusL(2),tangRadiusL(3),'k');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(AVO(1),AVO(2),AVO(3),...
                            LT(1),LT(2),LT(3),'k');
                q.AutoScaleFactor = 1;
                % ASSEMBLE VELOCITY OBSTACLE CONES
                nodes = 20;
                [Cone] = obj.vectorCone(AVO,Dvo,LTP,nodes);
            end
        end
        % VELOCITY OBSTACLES METHOD
        function [escapeVelocities] = getEscapeVelocities(obj,potentialVelocities,VOlist)
            % This function takes a matrix of potential velocities and
            % compares them against the agents Velocity Obstacle set. A
            % matrix of escape velocities is then produced.
            % INPUTS:
            % potentialVelocities - The feasible velocity matrix
            % VOlist              - The velocity obstacle list
            % OUTPUTS:
            % escapeVelocities    - The available escape velocity matrix
            
            escapeVelocities = [];   % Prepare output container
            
            % FOR EACH VIABLE VELOCITY, CHECK AGAINST CONES
            for i = 1:size(potentialVelocities,2)
                % Nominate a 3D velocity coordinate
                candidatePoint = potentialVelocities(:,i);
                % Check points against the VO set
                VOnumber = 1; flag = 0;
                while flag ~= 1 && VOnumber <= size(VOlist,2)
                    % Calculate axis properties
                    mod_coneAxis = sqrt(sum(VOlist(VOnumber).axis.^2));
                    unit_coneAxis = VOlist(VOnumber).axis/mod_coneAxis;
                    
                    % DETERMINE WHETHER POINT IS WITHIN VO
                    [flag] = obj.isInCone(VOlist(VOnumber).apex,unit_coneAxis,...
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
        
        % GET THE AGENT-SPECIFIC ENVIRONMENTAL UPDATE
        function [obj,obstacleSet,waypointSet] = getAgentUpdate(obj,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            % INPUT HANDLING
            if isempty(observedObjects)
                obstacleSet = [];
                waypointSet = [];
                return
            end
            % EXTRACT PARAMETERS TO BE STORED IN MEMORY
            for item = 1:length(observedObjects)
                % THE AGENT CAN ALLOCATE A PRIORITY TO OBJECTS WITHOUT ONE
                if ~isnan(observedObjects(item).priority)
                    objectPriority = observedObjects(item).priority;       % Assign priority to memory
                else
                    objectPriority = 1/observedObjects(item).range; 
                end
                % MEMORY ITEM CONSTRUCTOR
                memoryItem = struct('name',observedObjects(item).name,...
                                    'objectID',observedObjects(item).objectID,...
                                    'type',observedObjects(item).type,...
                                    'size',observedObjects(item).size,...
                                    'state',[observedObjects(item).position;observedObjects(item).velocity],...
                                    'inclinationAngle',observedObjects(item).inclinationAngle,...
                                    'azimuthAngle',observedObjects(item).azimuthAngle,...
                                    'priority',objectPriority); 
                % LOG TO MEMORY STRUCTURE
                obj = obj.updateAgentKnowledge(memoryItem);
            end
                
            % SORT OBJECT SET BY PRIORITY
            [~,ind] = sort([obj.memory.priority],2,'descend');             % Ordered indices of the object IDs
            obj.memory = obj.memory(ind);
            
            % DECERN OBSTACLES FROM WAYPOINTS
            % obj.memory now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            waypointIndices = [obj.memory.type] == simulation_objectType.waypoint; % Get the indices of the 'waypoint' types
            obstacleIndices = ~waypointIndices;
            waypointSet = obj.memory(waypointIndices);
            obstacleSet = obj.memory(obstacleIndices); 
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.waypointUpdater(waypointSet);                    % Update the target waypoint and heading vector
        end
    end
    %% STATIC MATHEMATIC FUNCTIONS
    methods (Static)
        % CONSTRUCT CONE FROM VECTOR DEFINTIONS
        function [Cone] = vectorCone(pointA,pointB,radialPoint,nodes)
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
            if ~exist('pointA','var')
                pointA = [0;0;0];
            end
            
            % MANUAL INPUTS
            coneColour = 'r';
            coneEdgeColour = 'r';
            coneAlpha = 0.1;
            coneClosed = 0;
            coneLines = 0;
            
            % DEFINE VECTOR PROBLEM PARAMETERS
            Doi = pointB - pointA;
            mDoi = sqrt(sum((Doi).^2)); % Axis vector properties
            tangent = radialPoint - pointA;
            mod_tangent = sqrt(sum((tangent).^2));       % Tangental vector properties
            
            % DEFINE THE TANGENT-AXIS PROJECTION
            trueAB = (dot(tangent,Doi)/mDoi^2)*Doi;
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
            angle_X1X2 = acos(dot(unit_Vx,Doi)/(norm(unit_Vx)*mDoi))*180/pi;
            
            % Get rotation axis
            axis_rot = cross([1 0 0],Doi);
            
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
        % COMPARE POINT TO CONE GEOMETRY
        function [flag] = isInCone(apex,unit_coneAxis,mod_coneAxis,coneAngle,point)
            % INPUTS:
            % Avo
            % mod_VOaxis
            % unit_VOaxis
            % alpha      - Cone angle
            % probePoint - Test point
            
            % Compare the angle made between the probeVector and the
            % unitVO axis. Is less than alpha? else not in cone.
            
            probeVector = point - apex;           % Define validation vector
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
        % DRAW SPHERE
        function [X,Y,Z] = getSphere(position,radius)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius
            [X,Y,Z] = sphere(40);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
        end
        % UNIT VECTOR
        function [unitVector] = unit(inputVector)
            % This function returns the unit vector of a vector
            unitVector = inputVector/sqrt(sum(inputVector.^2));
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]