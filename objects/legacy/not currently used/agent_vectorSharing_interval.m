%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing_interval.m) %%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_vectorSharing_interval < agent_vectorSharing
    %   AGENT PROPERTIES
    properties                  
        % INTERVAL UNCERTAINTY 
        confidenceAssumption = 3; % Standard deviations
    end
    %   CLASS METHODS
    methods
        %% CONSTRUCTOR
        function obj = agent_vectorSharing_interval(varargin)
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_vectorSharing(varargin);                             % Get super class 'agent'
            % CONTROL PARAMETERS
            obj.priorError = zeros(4,1);
%             obj.linearVelocityLimits = [30;30;1];                          % Limits on the agents velocity 
%             obj.linearAccelerationLimits = [1;1;1];                        % Limits on the agents acceleration
%             obj.angularVelocityLimits = [pi/6;pi/6;pi/6];
%             obj.angularAccelerationLimits = [0.1;0.1;0.1]; % rad/s^2
%             obj.angularAccelerationLimits = [0.05;0.1;0.1]; % rad/s^2
%             obj.angularVelocityLimits = [inf;inf;inf];
%             obj.angularAccelerationLimits = [0.05;0.1;0.1]; % rad/s^2
%             obj.angularAccelerationLimits = [inf;inf;inf]; % rad/s^2
            
            % SENSOR UNCERTAINTY
            obj.positionStandardDeviation = 0.5;          % Accurate to within 0.5m
            obj.velocityStandardDeviation = 0.1;          % Accurate to within 0.1m/s
            obj.rangeFinderStandardDeviation = 0.1;      % Accurate to within 0.2m
            obj.cameraStandardDeviation = 5.208E-5;       % 1 pixel in a 1080p image

            % DEFINE VIRTUAL SENSOR PROPERTIES            
            obj.VIRTUAL.size = 2.5;
            obj.VIRTUAL.detectionRange = obj.sensorRange; % Assign the range attribute to the SIM VIRTUAL attribute
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
            [obj,obstacleSet,~] = obj.getAgentUpdate(observationSet,dt);
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredSpeed = 5;
            if ~isempty(obj.targetWaypoint)  
               % GET THE WAYPOINT HEADING VECTOR
               desiredVelocity = obj.targetWaypoint.state(1:3)/norm(obj.targetWaypoint.state(1:3))*desiredSpeed;
            else
               desiredVelocity = [1;0;0]*desiredSpeed;
            end      
            
            % \\\\\\\\\\\\\\\\\\ OBSTACLE AVOIDANCE \\\\\\\\\\\\\\\\\\\\\\\
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(obstacleSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(desiredVelocity,obstacleSet,visualiseProblem);
            end
            algorithm_dt = toc(algorithm_start); 
           
            % \\\\\\\\\\\\\\\\ PARSE CONTROLLER INPUTS \\\\\\\\\\\\\\\\\\\\
            desiredVelocity = mid(desiredVelocity);
            targetSpeed = norm(desiredVelocity);
            targetHeading = desiredVelocity/targetSpeed;                            
            
            % \\\\\\\\\\\\\\\\\\\ AGENT CONTROLLER \\\\\\\\\\\\\\\\\\\\\\\\
            useController = 1;
            if useController
                [d_heading,d_speed,obj] = obj.trajectoryController(targetHeading,targetSpeed);  
            else
                [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);
            end
            newSpeed = (norm(obj.localState(7:9)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = d_heading;                                        % Heading change relative to straight ahead  
            
            % \\\\\\\\\\\\\\ APPLY KINEMATIC LIMITS \\\\\\\\\\\\\\\\\\\\\\
            [newHeading,newSpeed] = obj.kinematicContraints(dt,newHeading,newSpeed);
            newVelocity = [1;0;0]*newSpeed;                                % Convert to local velocity vector
                     
            % \\\\\\\\\\\\\\\ GET THE NEW STATE VECTOR \\\\\\\\\\\\\\\\\\\\
            newState = obj.stateDynamics_velocities(dt,newVelocity,newHeading);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                newState = obj.freezeAgent();
            end  
            % \\\\\\\\\\\ UPDATE THE CLASS GLOBAL PROPERTIES \\\\\\\\\\\\\\
            obj = obj.updateGlobalProperties(dt,newState);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:length(obj.priorError),TIME.currentStep) = [newSpeed;newState(4:6)];       % Record the control inputs
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
        end
        
        % GET THE INTERVAL SHARING VELOCITY CORRECTION
        function [avoidanceVector] = getAvoidanceCorrection(obj,desiredVelocity,obstacleSet,visualiseProblem)
            % This function is the highest level of the obstacle avoidance
            % algorithm. It accepts the agents current desired trajectory
            % and obstacles the agent is aware of (currently in memory) and
            % computes an alternative trajectory based on the given algorithm.
             
            avoidanceHorizon = 10;
            
            % UPDATE (MEASURE) OWN TRAJECTORY
            [~,absVelocityA,zeroPositionA,zeroVelocityA] = obj.getAgentMeasurements(); % Get the uncertain position and uncertainties
            
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            avoidanceVelocitySet = [];
            for item = 1:length(obstacleSet)
                % GET THE OBSTACLES INFORMATION
                radiusIntvalB   = obstacleSet(item).size;                  % Get object size from memory
                positionIntvalB = obstacleSet(item).state(1:3);            % Get current position from memory
                velocityIntvalB = obstacleSet(item).state(4:6);            % Get current velocity from memory
                
                % CHECK WHETHER OBSTACLE PROXIMITY IS VALID FOR AVOIDANCE
                if inf(norm(positionIntvalB)) > avoidanceHorizon
                    continue                                               % Continue if obstacle is not to be considered
                end
                                            
                % GET COLLISION LIKELYHOOD
                [tau] = obj.validateCollision(positionIntvalB,velocityIntvalB);
                if any(isnan(velocityIntvalB)) || inf(tau) < 0             % Only abort if all of tau is negative
                    continue                                               % Avoidance requires a velocity estimate
                end
                
                % COMPUTE THE RELATIVE POSITION AND VELOCITY WITH ADDED
                % UNCERTAINTY IN THE CURRENT AGENTS POSITION AND VELOCITY
                positionEstimateB = positionIntvalB - zeroPositionA;       % S = S_b - S_a
                velocityEstimateB = velocityIntvalB - zeroVelocityA;       % V = V_b - V_a
                
                % COMPUTE THE VECTOR SHARING PROBLEM
                [Voptimal] = obj.defineProblem(absVelocityA,positionEstimateB,...
                                               velocityEstimateB,radiusIntvalB,visualiseProblem);
                
                Voptimal                            
                % CONCATINATE THE OPTIMAL SETS
                if ~isnan(Voptimal)
                    avoidanceVelocitySet = horzcat(avoidanceVelocitySet,Voptimal);
                end
            end

            % \\\\\\\\\\\\\\ COMPUTE OBSTACLE INTERSECTIONS \\\\\\\\\\\\\\
            if ~isempty(avoidanceVelocitySet)
                % CALCULATE THE GLOBAL OPTIMAL HEADING
                avoidanceVector = obj.evaluateAvoidanceIntersections(desiredVelocity,avoidanceVelocitySet);
            else
                avoidanceVector = desiredVelocity;
            end
            % CONVERT TO AN INACTABLE VELOCITY 
            avoidanceVector = mid(avoidanceVector);
        end
        % DEFINE THE INTERVAL GEOMETRIC AVOIDANCE PROBLEM %%%%%%%%%%%%%%%%%
        function [U_a] = defineProblem(obj,V_a,S_b,c_b,obstacleSize,visualiseProblem)
            % This function calculates the avoidance vectors based on the
            % principle of vector sharing. 
            % S_a - Estimate of the local position
            % V_a - Estimate of the absolute velocity
            % S_b - Estimate of the relative obstacle position
            % c_b - Estimate of the relative obstacle velocity
            % obstacleSize - Estimate of the obstacle size
            
            % IN THE EVENT THE INPUT VELOCITY IS ZERO     
            
            % DEFINE THE INITIAL PROBLEM PARAMETERS
            r = S_b;                                  % The relative position vector = relative seperation
            V_b = V_a + c_b;                          % Define absolute velocity of B 
            c_b_unit = obj.intervalUnit(c_b);         % The unit relative velocity of B
           
            % DEFINE THE MISS VECTOR           
            temp = obj.intervalCross(r,c_b_unit);
            r_m  = obj.intervalCross(c_b_unit,temp);                                    % The near-miss vector
            
            tol = 1E-5;
            if norm(mid(r_m)) < tol 
                temp = obj.intervalCross(r,c_b_unit);
                r_m  = obj.intervalCross(c_b_unit,temp); 
                r_m = r_m + rand(3,1);
            end

            % GET THE NORM OF THE MISS INTERVAL
            rNorm = obj.intervalNorm(r_m);
            
            % DEFINE THE TIME TO CLOSEST APPROACH(+ve converging)
            tau = -(dot(r,c_b)/dot(c_b,c_b));

            % DEFINE THE REST REGION
            r_safe = obj.VIRTUAL.size + obj.obstacleSafetyFactor*(obstacleSize); % The specified minimum safe seperation limit (m)       
            r_res  = r_safe - rNorm;
            
            % RESOLUTION RATION
            resRatio = r_res/rNorm;
            % DEFINE THE VECTOR SHARING TERMS
%             r_vsa = (norm(V_b)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(-r_m)
%             r_vsb = (norm(V_a)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(r_m)
            AShareRatio = obj.intervalNorm(V_b)/(obj.intervalNorm(V_a) + obj.intervalNorm(V_b));  
            r_vsa = AShareRatio*resRatio*-r_m;
                        
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (V_a*tau + r_vsa);
        end
        
        % //////////////////// SENSOR FUNCTIONS ///////////////////////////
        % AGENT-SPECIFIC UPDATE FUNCTION (USING SPHERICAL DATA & SENSOR MODEL)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,observationSet,dt)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observationSet - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            for item = 1:length(observationSet)
                % DEFINE THE SENSED VARIABLES
                % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
                objectName = observationSet(item).name;                    % The objects name
                objectID   = observationSet(item).objectID;                % The objects ID
                objectType = observationSet(item).type;                    % The objects type (obstacle or waypoint)
                objectPriority = observationSet(item).priority;            % The objects priority
                % ADDITIONAL MEASUREMENTS
                [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.getSensorMeasurments(observationSet(item));
                [measuredPosition] = obj.cartesianFromSpherical(measuredRange,measuredAzimuth,measuredElevation);    % Calculate new relative position
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                objectSize = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % GET THE AGENTS KNOWLEDGE OF THIS OBSTACLE
                [priorState] = obj.getLastStateFromMemory(objectID);       % Get prior obstacle knowledge
                if isempty(priorState)
                    % FIRST SIGHT OF OBSTACLE
                    objectState = [measuredPosition;NaN(3,1)];             % If no prior knowledge
                else
                    priorState = mid(priorState);                          % Assume the previous measurment was correct
                    % CALCULATE THE STATE UPDATE
                    [objectState] = obj.linearStateEstimation(priorState,measuredPosition,dt);
                end
                
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(objectPriority)
                    objectPriority = observationSet(item).priority;           % Assign priority to memory
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
            if isintval([obj.memory.priority])
                [~,ind] = sort(mid([obj.memory.priority]),2,'descend');             % Ordered indices of the object IDs
                obj.memory = obj.memory(ind);
            else
                [~,ind] = sort([obj.memory.priority],2,'descend');             % Ordered indices of the object IDs
                obj.memory = obj.memory(ind);
            end
            % DECERN OBJECT TYPES
            % obj.memory now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            agentIndices = [obj.memory.type] == simulation_objectType.agent;
            waypointIndices = [obj.memory.type] == simulation_objectType.waypoint;  
            obstacleIndices = [obj.memory.type] == simulation_objectType.obstacle; % Differenciate between the different object types
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet = obj.memory(agentIndices);
            waypointSet = obj.memory(waypointIndices);
            obstacleSet = obj.memory(obstacleIndices); 
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.waypointUpdater(waypointSet);  % Update the target waypoint and heading vector
        end
        % AGENT-SPECIFIC SELF MEASUREMENTS
        function [absPositionA,absVelocityA,zeroPositionA,zeroVelocityA] = getAgentMeasurements(obj)
            % GET THE STATE UNCERTAINTIES
            positionAUncertainty = obj.positionStandardDeviation*randn(3,1);                    % Zero estimate of position
            velocityAUncertainty = obj.velocityStandardDeviation*randn(3,1);                    % Zero estimate of velocity
            % GET THE LOCAL RELATIVE VARIABLES
            zeroPositionA = midrad(positionAUncertainty,3*obj.positionStandardDeviation);
            zeroVelocityA = midrad(velocityAUncertainty,3*obj.velocityStandardDeviation); 
            % GET THE LOCAL ABSOLUTE VARIABLES
            absPositionA = zeroPositionA + obj.localState(1:3); % Get the absolute position measurement
            absVelocityA = zeroVelocityA + obj.localState(4:6); % Get the absolute velocity measurement          
        end
        % AGENT-SPECIFIC SENSOR MODEL
        function [rangeBox,azimuthBox,elevationBox,alphaBox] = getSensorMeasurments(obj,obstacleData) 
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
            measuredRange     = obstacleData.range + obj.rangeFinderStandardDeviation*randn(1);     
            measuredAzimuth   = obstacleData.azimuthAngle + obj.cameraStandardDeviation*randn(1);
            measuredElevation = obstacleData.inclinationAngle + obj.cameraStandardDeviation*randn(1);
            % OBSERVED ANGULAR WIDTH
            measuredAlpha = obstacleData.angularWidth + obj.cameraStandardDeviation*randn(1);      % Uncertainty in the angular measurement
            
            %% DEFINE THE MEASUREMENT INTERVALS
            rangeBox     = midrad(measuredRange,obj.confidenceAssumption*obj.rangeFinderStandardDeviation);
            azimuthBox   = midrad(measuredAzimuth,obj.confidenceAssumption*obj.cameraStandardDeviation);
            elevationBox = midrad(measuredElevation,obj.confidenceAssumption*obj.cameraStandardDeviation);
            alphaBox     = midrad(measuredAlpha,obj.confidenceAssumption*obj.cameraStandardDeviation);
        end
    end
%   STATIC METHODS    
    methods (Static) % Don't require an input object to work        
        % GET THE PITCH AND HEADING ERROR FROM THE INTERVAL INPUTS
        function [theta,lambda] = getControlInputs(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS
            % V - The current unity velocity vector
            % U - The desired heading vector
            
            V = mid(V);
            U = mid(U);
            
            %% GET THE LINE OF SIGHT ANGLE             
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = diag([1 1 0])*V;
            Uh = diag([1 1 0])*U;                                          % Reject the vertical elements
            % DEFINE LOS ANGLE                        
            rotationAxis = cross(Vh,Uh);
            rotationSign = sign(mid(rotationAxis(3)));        
            % CALCULATE THE RELATIVE HEADING ANGLE
            lambda = rotationSign*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            if isnan(lambda)
                lambda = 0;
            end
            % GET THE ELEVATION ANGLE
            if isintval(U)
                theta = atan(U(3)/norm(Uh));
            else
                theta = atan2(U(3),norm(Uh));
            end
            if isnan(theta)
                theta = 0;
            end
        end
        % CALCULATE THE GLOBAL OPTIMAL INTERVAL
        function [optimalHeading] = evaluateAvoidanceIntersections(waypointHeading,avoidanceVelocities)
            % This function calculates the target heading from a set of
            % avoidance velocities. The intersection returns [nan,nan,nan]
            % in the event of a failed intersection. The last successful
            % intersection should ideally be returned in this case.
            % INPUT:
            % avoidanceVelocities - The list of avoidance intervals
            % OUTPUT:
            % targetHeading       - The resulting global (singular) heading
            
            % INPUT HANDLING
            if isempty(avoidanceVelocities)
               optimalHeading = mid(waypointHeading);
               return
            end
            
            
            % MOVE THROUGH THE VELOCITY SET
            viableSet = avoidanceVelocities(:,1);                          % Initialise the global intersection with the highest priority           
            for item = 1:size(avoidanceVelocities,2)
                % GET THE NEW AVOIDANCE INTERVAL
                avoidanceInterval = avoidanceVelocities(:,item);           
                
                % INTERSECT THE CURRENT INTERVAL AGAINST THE GLOBAL SET
                globalIntersection = intersect(viableSet,avoidanceInterval); 
                if ~any(isnan(globalIntersection))
                    % THE NEW VIABLE SET IS THE INTERSECTION
                    viableSet = globalIntersection;
                end
            end
            % CALCULATE THE UNIT HEADING VECTOR FROM THE GLOBAL INTERVAL
%             optimalHeading = avoidanceVelocities(:,1); 
            optimalHeading = mid(globalIntersection);
            
            % ATTEMPT INTERSECTION WITH WAYPOINT HEADING
%             optimalHeading = intersect(globalIntersection,waypointHeading)
%             if any(isnan(optimalHeading))
%                 % IF FAILED, ASSUME AVOIDANCE HEADING
%                 optimalHeading = globalIntersection;
%             end
            % CALCULATE THE UNIT HEADING VECTOR FROM THE GLOBAL INTERVAL
%             optimalHeading = mid(optimalHeading);
        end
        % VALIDATE THE OBSTACLE
        function [tau] = validateCollision(S_b,c_b)
            % CONFIRM WE KNOW THE HEADING OF THE OBSTACLE
            if any(isnan(c_b))
                tau = -inf;
                return
            end                 
            % Define the time to closest approach (+ve converging)
            A = dot(S_b,c_b)
            B = dot(c_b,c_b)
            tau = -(A/B)
        end
        
        % MATH TOOLS
        % CONVERT CARTESIAN TO SPHERICAL
        function [range,azimuth,elevation] = sphericalFromCartesian(positionVector)
            range = norm(positionVector);                                            % The range
%             azimuth = atan2(positionVector(2),positionVector(1));                    % The angle made in the azimuth (bearing)
%             elevation = atan2(positionVector(3),sqrt(positionVector(1).^2 + positionVector(2).^2));  % The elevation (vertical bearing)
            azimuth = atan(positionVector(2)/positionVector(1))
            elevation = atan(positionVector(3)/sqrt(positionVector(1).^2 + positionVector(2).^2))
        end
        % INTERVAL CROSS PRODUCT
        function [crossInt] = intervalCross(intA,intB)
            crossInt = [intA(2).*intB(3)- intA(3).*intB(2); 
                        intA(3).*intB(1)- intA(1).*intB(3); 
                        intA(1).*intB(2)- intA(2).*intB(1)]; 
        end
        % INTERVAL UNIT VECTOR
        function [unitInt] = intervalUnit(intVector,zeroOmitted)
            % INPUT HANDLING
            if ~exist('zeroOmitted','var')
                zeroOmitted = 0;
            end
            % CHECK IF NOT INTERVAL TYPE
            if ~isintval(intVector)
                % If the input value is not an interval
                unitInt = intVector/norm(intVector);
                return
            end
            % IF ZERO NOT IN RANGE
            if zeroOmitted
                for i = 1:size(intVector)
                    s = 0; % Collect the square terms
                    for j = 1:size(intVector)
                        if i == j
                            continue
                        end
                        s = s + intVector(j)^2;
                    end
                    x = intVector(i);
                    
                    % COMPUTE THE INTERVAL NORM WHERE 0 IS NOT INCLUDED
                    n(i) = sign(inf(x))/sqrt(1+(s/x^2));
                end
            else
                % ZERO IS PART OF THE RANGE ( GENERALISED )
                for i = 1:size(intVector,1)
                    s = 0; % Collect the square terms
                    for j = 1:size(intVector,1)
                        if i == j
                            continue
                        end
                        s = s + sqr(intVector(j));
                    end
                    x = intVector(i);
                    
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
            unitInt = n;
            % CONFIRM FORMATTING OF THE VECTOR
            unitInt = reshape(unitInt,size(intVector));
        end
        % INTERVAL NORM VECTOR
        function [normInt] = intervalNorm(intVector)
            % This function calculates the norm of an interval vector
            % CHECK IF NOT INTERVAL TYPE
            if ~isintval(intVector)
                % If the input value is not an interval
                normInt = norm(intVector);
                return
            end
%             sqrVector = intVector.^2;
%             sumVector = sum(sqrVector);
%             normInt = sqrt(sumVector);
            
%             offset = 0.2E-2;
%             if inf(normInt) == 0
%                 normInt = infsup(offset ,sup(normInt));
%             end

            infVal = inf(intVector).^2;
            supVal = sup(intVector).^2;
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
            normInt = infsup(sqrtMatrix(1),sqrtMatrix(2));     
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]