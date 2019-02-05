%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_vectorSharing < agent
    
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
%   AGENT PROPERTIES 
    properties     
        
        % ALGORTHIM PROPERTIES
        obstacleSafetyFactor = 4;         % Modify the apparent size of the obstacle
        % SENSOR ACCURACY PARAMETERS
        % SELF-SENSING
        positionStandardDeviation;          % Accurate to within 0.5m
        velocityStandardDeviation;          % Accurate to within 0.1m/s
        % VISUAL SENSING
        rangeFinderStandardDeviation;       % Accurate to within 0.2m
        cameraStandardDeviation;            % 1 pixel in a 1080p image  
        
        % PID TERMS
        priorError = zeros(3,1);
        
        % Giff Toggle
        giffOn = 0;
    end
%   CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_vectorSharing(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'
            % PARAMETERISE KINEMATIC CONSTRAINTS
            obj.linearVelocityLimits = [30;30;1];                          % Limits on the agents velocity 
            obj.linearAccelerationLimits = [1;1;1];                        % Limits on the agents acceleration
            obj.angularVelocityLimits = [pi/3;pi/3;pi/3];
            obj.angularAccelerationLimits = [0.05;0.2;0.2]; % rad/s^2
            % SENSOR ATTRIBUTES
            obj.positionStandardDeviation = 0;                             % Accurate to within 0.5m
            obj.velocityStandardDeviation = 0;                             % Accurate to within 0.1m/s
            obj.rangeFinderStandardDeviation = 0;                          % Accurate to within 0.2m
            obj.cameraStandardDeviation = 0;                               % 1 pixel in a 1080p image  
            obj.sensorRange = 200;                                         % Sensor range parameter
            % DEFINE VIRTUAL PROPERTIES
            obj.VIRTUAL.detectionRange = obj.sensorRange;                  % Assign the range attribute to the SIM VIRTUAL attribute     
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        %% AGENT MAIN CYCLE 
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
                   
            % CHECK FOR NEW INFORMATION UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\\\\
            observationSet = varargin{1}; % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,~] = obj.getAgentUpdate(observationSet,dt);

            % GET WAYPOINT (TARGET) HEADING VECTOR \\\\\\\\\\\\\\\\\\\\\\\\
            if ~isempty(obj.targetWaypoint)  
               % GET THE WAYPOINT HEADING VECTOR
               targetHeading = obj.targetWaypoint.state(1:3)/norm(obj.targetWaypoint.state(1:3));
            else
               targetHeading = [1;0;0];
            end

            % \\\\\\\\\\\ UPDATE (MEASURE) OWN TRAJECTORY \\\\\\\\\\\\\\\\\                               
            [~,velocityEstimateA,zeroPositionA,zeroVelocityA] = obj.getAgentMeasurements();
 
            % \\\\\\\\\\\ OBSTACLE AVOIDANCE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            algorithm_start = tic; algorithm_indicator = 0;    
            avoidanceVelocitySet = []; visualiseAgent = 1; visualiseProblem = 0;  
            if ~isempty(obstacleSet)
                % MOVE THROUGH THE PRIORITISED OBSTACLE SET
                for item = 1:length(obstacleSet)
                    % GET THE OBSTACLES INFORMATION
                    sizeEstimateB     = obstacleSet(item).size;       % Get object size from memory
                    positionEstimateB = obstacleSet(item).state(1:3); % Get current position from memory
                    velocityEstimateB = obstacleSet(item).state(4:6); % Get current velocity from memory
                    % GET COLLISION LIKELYHOOD
                    [tau] = obj.validateCollision(positionEstimateB,velocityEstimateB);
                    if any(isnan(velocityEstimateB)) || tau < 0
                       continue                                       % Avoidance requires a velocity estimate
                    end
                    % INDICATE ALGORITHM WAS RAN
                    algorithm_indicator = 1;
                    % COMPUTE THE RELATIVE POSITION AND VELOCITY WITH ADDED
                    % UNCERTAINTY IN THE CURRENT AGENTS POSITION AND VELOCITY
                    positionEstimateB = positionEstimateB - zeroPositionA;   % S = S_b - S_a
                    velocityEstimateB = velocityEstimateB - zeroVelocityA;   % V = V_b - V_a
                    % COMPUTE THE VECTOR SHARING PROBLEM
                    [Voptimal,~] = obj.defineProblem(velocityEstimateA,positionEstimateB,...
                                                     velocityEstimateB,sizeEstimateB,...
                                                     visualiseAgent,visualiseProblem);
                    % CONCATINATE THE OPTIMAL SETS
                    avoidanceVelocitySet = horzcat(avoidanceVelocitySet,Voptimal);
                end
            end
            % IF THERE ARE AVOIDANCE TRAJECTORIES, SELECT HIGHEST PRIORITY
            if ~isempty(avoidanceVelocitySet) 
                % SELECT THE HIGHEST PRIORITY OBSTACLE
                optimalHeading = avoidanceVelocitySet(:,1);
                targetHeading = optimalHeading/norm(optimalHeading);
            end
            algorithm_dt = toc(algorithm_start);
            
            % \\\\\\\\\\\\\\ CALCULATE CONTROL INPUTS \\\\\\\\\\\\\\\\\\\\\
            [theta,lambda] = obj.getControlInputs([1;0;0],targetHeading);
            controlError = [0;theta;-lambda];           % -ve in the NED control frame
            % BUILD THE CONTROL INPUTS
            Kp = 0.3; Kd = 10;    % 0.2 , 10
            ctrlAcceleration = Kp*eye(3)*controlError + Kd*eye(3)*(controlError-obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = controlError;

            % \\\\\\\\\\\\\\ APPLY KINEMATIC LIMITS \\\\\\\\\\\\\\\\\\\\\\\
            linearAcceleration = zeros(3,1);
            % BOUND ANGULAR ACCELERATION VALUES
            [ctrlAcceleration(1)] = obj.boundValue(ctrlAcceleration(1),-obj.angularAccelerationLimits(1),obj.angularAccelerationLimits(1));
            [ctrlAcceleration(2)] = obj.boundValue(ctrlAcceleration(2),-obj.angularAccelerationLimits(2),obj.angularAccelerationLimits(2));
            [ctrlAcceleration(3)] = obj.boundValue(ctrlAcceleration(3),-obj.angularAccelerationLimits(3),obj.angularAccelerationLimits(3));
            % \\\\\\\\\\\\ GET THE NEW STATE VECTOR \\\\\\\\\\\\\\\\\\\\\\\
            agentStateUpdate = obj.stateDynamics(dt,linearAcceleration,ctrlAcceleration);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                agentStateUpdate = obj.freezeAgent();
            end 
            % \\\\\\\\\\\\ UPDATE THE CLASS GLOBAL PROPERTIES \\\\\\\\\\\\\
            obj = obj.updateGlobalProperties(dt,agentStateUpdate);
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.controlVector(1:3,TIME.currentStep) = ctrlAcceleration;      % Record the control inputs 
        end
        
        % DEFINE THE GEOMETRIC AVOIDANCE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [U_a,U_b,unitAcceleration] = defineProblem(obj,V_a,S_b,c_b,obstacleSize,visualiseAgent,visualiseProblem)
            % This function calculates the avoidance vectors based on the
            % principle of vector sharing. 
            % S_a - Estimate of the local position
            % V_a - Estimate of the local velocity
            % S_b - Estimate of the relative obstacle position
            % c_b - Estimate of the relative obstacle velocity
            % obstacleSize - Estimate of the obstacle size

%             % IN THE EVENT THE INPUT VELOCITY IS ZERO
%             if norm(c_b) == 0 
%                 % INSUFFICIENT VELOCITY DATA
%                 U_a = zeros(3,1);
%                 U_b = zeros(3,1);
%                 return
%             end

            % TO PREVENT A ZERO MISS VECTOR
            if ~any(cross(S_b,c_b))
                S_b = S_b + randn(3,1);
            end
            
            % DEFINE THE INITIAL PROBLEM PARAMETERS
            r = S_b;                                  % The relative position vector = relative seperation
            V_b = V_a + c_b;                          % Define absolute velocity of B 
            c_b_unit = c_b/norm(c_b);                 % Unit relative velocity of V_b  
            r_m = cross(c_b_unit,cross(r,c_b_unit));  % The near-miss vector 
            % DEFINE THE TIME TO COLLISION (+ve converging)
            tau = -(dot(S_b,c_b)/dot(c_b,c_b));                                       
            % DEFINE THE REST REGION
            r_safe = obj.VIRTUAL.size + obj.obstacleSafetyFactor*(obstacleSize);   % The specified minimum safe seperation limit (m)
            r_res = r_safe - norm(r_m);                                            % Define the rest region
            % DEFINE THE VECTOR SHARING TERMS
            r_vsa = (norm(V_b)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(-r_m);
            r_vsb = (norm(V_a)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(r_m);
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (V_a*tau + r_vsa);
%             Unorm_a = U_a/(norm(V_a*tau + r_vsa));
            U_b = (V_b*tau + r_vsb);
%             Unorm_b = U_b/(norm(V_b*tau + r_vsb));
            
            % DEFINE THE OPTIMAL ACCELERATION VECTOR
            % Based on the hamiltonian, H, the vector that minimises the
            % cost expression J, is an ahat that is parallel to -r_m:
            % J = -0.5*norm(r_m)^2;
            % H = dot(-r_m,c_b) + (a_mod*tau)*dot(-r_m,unit_acc);
            
            % We can therefore write the optimal acceleration direction:
            unitAcceleration = -r_m/norm(r_m);

            % PLOT THE SCENARIO IF REQUESTED
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                % A & B's current trajectory at tau
                Vtau_a = V_a*tau;
                Vtau_b = V_b*tau;
                
                % The vector sharing parameters
                temp = Vtau_b + S_b;
                
                hold on; grid on;
                axis equal;
                
                % DRAW THE INITIAL PROBLEM
                q = quiver3(0,0,0,V_a(1),V_a(2),V_a(3),'r'); q.AutoScaleFactor = 1;
                q = quiver3(S_b(1),S_b(2),S_b(3),V_b(1),V_b(2),V_b(3),'b'); q.AutoScaleFactor = 1;
                q = quiver3(S_b(1),S_b(2),S_b(3),c_b(1),c_b(2),c_b(3),'m'); q.AutoScaleFactor = 1; 

                % FIGURE 3
                % A & B's current trajectory at tau
                q = quiver3(0,0,0,Vtau_a(1),Vtau_a(2),Vtau_a(3),'k'); q.AutoScaleFactor = 1;          % A's current trajectory at tau
                q = quiver3(r(1),r(2),r(3),Vtau_b(1),Vtau_b(2),Vtau_b(3),'k'); q.AutoScaleFactor = 1; % B's current trajectory at tau

                % Near-miss vector @ tau
                q = quiver3(Vtau_a(1),Vtau_a(2),Vtau_a(3),r_m(1),r_m(2),r_m(3),'k'); q.AutoScaleFactor = 1;
                % The vector sharing parameters
                q = quiver3(Vtau_a(1),Vtau_a(2),Vtau_a(3),r_vsa(1),r_vsa(2),r_vsa(3),'r'); q.AutoScaleFactor = 1;
                q = quiver3(temp(1),temp(2),temp(3),r_vsb(1),r_vsb(2),r_vsb(3),'b'); q.AutoScaleFactor = 1;
                % The desired correction vectors
                q = quiver3(0,0,0,U_a(1),U_a(2),U_a(3),'g'); q.AutoScaleFactor = 1;
                q = quiver3(r(1),r(2),r(3),U_b(1),U_b(2),U_b(3),'g'); q.AutoScaleFactor = 1;
                
                %% DEFINE THE RESOLUTION ZONE
%                 zoneCenter = Vtau_a + 0.5*r_m;
%                 zoneRadius = 0.5*(norm(r_vsa) + norm(r_vsb) + norm(r_m));  % Define the radius of the safe region
%                 % DRAW THE RESOLUTION ZONE
%                 scatter3(zoneCenter(1),zoneCenter(2),zoneCenter(3),'k');
%                 [X_zone,Y_zone,Z_zone] = obj.getSphere(zoneCenter,zoneRadius);
%                 sphZone = mesh(X_zone,Y_zone,Z_zone);
%                 set(sphZone,'facealpha',0.2,...
%                             'FaceColor','g',...
%                             'LineWidth',0.1,...
%                             'EdgeAlpha',0.2);
                
                %% DRAW THE POSITION ADN SIZE OF THE AGENTS
                % DRAW THE RESOLUTION ZONE
                [X_A,Y_A,Z_A] = obj.getSphere([0;0;0],obj.VIRTUAL.size);
                sphZone = mesh(X_A,Y_A,Z_A);
                set(sphZone,'facealpha',0.2,...
                            'FaceColor','g',...
                            'LineWidth',0.1,...
                            'EdgeAlpha',0.2);
                [X_A,Y_A,Z_A] = obj.getSphere(r,obstacleSize);
                sphZone = mesh(X_A,Y_A,Z_A);
                set(sphZone,'facealpha',0.2,...
                            'FaceColor','g',...
                            'LineWidth',0.1,...
                            'EdgeAlpha',0.2);
                %% DRAW THE ACCELERATION
                q = quiver3(0,0,0,unitAcceleration(1),unitAcceleration(2),unitAcceleration(3),'r'); q.AutoScaleFactor = 1;
            end
        end

        % AGENT-SPECIFIC UPDATE FUNCTION (USING SPHERICAL DATA & SENSOR MODEL)
        function [obj,obstacleSet,waypointSet] = getAgentUpdate(obj,observedObjects,dt)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            for item = 1:length(observedObjects)
                %% DEFINE THE SENSED VARIABLES
                % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
                objectName = observedObjects(item).name;                   % The objects name
                objectID = observedObjects(item).objectID;                 % The objects ID
                objectType = observedObjects(item).type;                   % The objects type (obstacle or waypoint)
                objectPriority = observedObjects(item).priority;           % The objects priority
                [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.getSensorMeasurments(observedObjects(item));
                [measuredPosition] = obj.cartesianFromSpherical(measuredRange,measuredAzimuth,measuredElevation);    % Calculate new relative position
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                objectSize = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % GET THE AGENTS KNOWLEDGE OF THIS OBSTACLE
                [priorState] = obj.getLastStateFromMemory(objectID);   % Get prior obstacle knowledge
                if isempty(priorState)
                    % FIRST SIGHT OF OBSTACLE
                    objectState = [measuredPosition;NaN(3,1)];         % If no prior knowledge
                else                   
                    % CALCULATE THE STATE UPDATE
                    [objectState] = obj.linearStateEstimation(priorState,measuredPosition,dt);
                end
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(objectPriority)
                    objectPriority = observedObjects(item).priority;           % Assign priority to memory
                else
                    objectPriority = 1/measuredRange; 
                end
                % MEMORY ITEM CONSTRUCTION
                memoryItem = struct('name',objectName,'objectID',objectID,...
                                    'type',objectType,'range',measuredRange,...
                                    'azimuth',measuredAzimuth,'elevation',measuredElevation,...
                                    'alpha',measuredAlpha,'size',objectSize,...
                                    'state',objectState,'priority',objectPriority); 
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.updateAgentKnowledge(memoryItem);
            end
            % SORT OBJECT SET BY PRIORITY
            priorityVector = [obj.memory.priority];
            [~,ind] = sort(priorityVector,2,'descend');             % Ordered indices of the object IDs
            obj.memory = obj.memory(ind);
            % DECERN OBSTACLES FROM WAYPOINTS
            % obj.memory now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            waypointIndices = strcmpi({obj.memory.type},'WAYPOINT');       % Get the indices of the 'waypoint' types
            obstacleIndices = ~waypointIndices;
            waypointSet = obj.memory(waypointIndices);
            obstacleSet = obj.memory(obstacleIndices);  
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.waypointUpdater(waypointSet);  % Update the target waypoint and heading vector
        end
        % AGENT SELF MEASUREMENTS
        function [absPositionA,absVelocityA,zeroPositionA,zeroVelocityA] = getAgentMeasurements(obj)
            % GET THE LOCAL RELATIVE VARIABLES
            zeroPositionA = zeros(3,1);
            zeroVelocityA = zeros(3,1); 
            % GET THE LOCAL ABSOLUTE VARIABLES
            absPositionA = zeroPositionA + obj.localState(1:3); % Get the absolute position measurement
            absVelocityA = zeroVelocityA + obj.localState(4:6); % Get the absolute velocity measurement          
        end
        % AGENT SENSOR MODEL
        function [range,azimuth,elevation,angularWidth] = getSensorMeasurments(obj,obstacleData) 
            % This function takes the simulation data and calculates the
            % spherical position and radius that would otherwise be sensed 
            % by the system.
            % INPUTS:
            % obstacleData - The object observation structure
            % OUTPUTS:
            % range        - The obstacles apparent range
            % azimuth      - The obstacles angular position in the XY plane
            % elevation    - The obstacles angular position in the XZ plane
            % angularWidth - The obstacles angular width in the azimuth

            % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
            range     = obstacleData.range + obj.rangeFinderStandardDeviation*randn(1);     
            azimuth   = obstacleData.azimuthAngle + obj.cameraStandardDeviation*randn(1);
            elevation = obstacleData.inclinationAngle + obj.cameraStandardDeviation*randn(1);
            % OBSERVED RADIUS
            angularWidth = obstacleData.angularWidth + obj.cameraStandardDeviation*randn(1);      % Uncertainty in the angular measurement
        end
    end
%   STATIC METHODS
    methods (Static)
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]