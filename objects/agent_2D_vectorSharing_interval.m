%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing_interval.m) %%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_2D_vectorSharing_interval < agent_vectorSharing_interval & agent_2D_vectorSharing
    %   AGENT PROPERTIES
    properties                  
        % INTERVAL UNCERTAINTY 
%         confidenceAssumption = 3; % Standard deviations
    end
    %   CLASS METHODS
    methods
        %% CONSTRUCTOR
        function obj = agent_2D_vectorSharing_interval(varargin)
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_vectorSharing(varargin);                             % Get super class 'agent'
            
%             % SENSOR UNCERTAINTY
%             obj.positionStandardDeviation = 0.5;          % Accurate to within 0.5m
%             obj.velocityStandardDeviation = 0.1;          % Accurate to within 0.1m/s
%             obj.rangeFinderStandardDeviation = 0.1;      % Accurate to within 0.2m
%             obj.cameraStandardDeviation = 5.208E-5;       % 1 pixel in a 1080p image

            % GET THE SENSOR DESCRIPTIONS
            obj = obj.getDefaultSensorParameters();     % For defaults
%             obj = obj.getCustomSensorParameters();    % For experiment

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        %% AGENT MAIN CYCLE
        function [obj] = main(obj,TIME,varargin)
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
            
            % /////////////////// DEFAULT BEHAVIOUR ///////////////////////
            desiredSpeed = obj.maxSpeed;
            desiredHeadingVector = [1;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(dt,varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
                desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
            end
 
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = obj.getAvoidanceCorrection(desiredVelocity,avoidanceSet,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % LIMIT THE SPEED
            if desiredSpeed > obj.maxSpeed
                desiredSpeed = obj.maxSpeed;
            end
            % TAKE THE CENTER OF THE INTERVAL FOR ACTUATION
            desiredVelocity = mid(desiredHeadingVector*desiredSpeed);
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [newState] = obj.stateDynamics_simple(dt,desiredVelocity,0);
            obj = obj.updateGlobalProperties_fixedFrame(dt,newState);  
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputNames = {'Vx (m/s)','Vy (m/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = newState(4:6);         % Record the control inputs
            
            % // DISPLAY CONFLICT RESOLUTION
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                obj = obj.getAnimationFrame(overHandle,TIME,'resolutionZone.gif');
                close(overHandle);
            end
        end
        
        % GET THE INTERVAL SHARING VELOCITY CORRECTION
        function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,knownObstacles,visualiseProblem)
            % This function is the highest level of the obstacle avoidance
            % algorithm. It accepts the agents current desired trajectory
            % and obstacles the agent is aware of (currently in memory) and
            % computes an alternative trajectory based on the given algorithm.

            % AGENT KNOWLEDGE
            [p_a,v_a,r_a] = obj.getAgentMeasurementIntervals(); % Its own position, velocity and radius
           
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            optimalSet = [];
            for item = 1:length(knownObstacles)
                % CHECK NEIGHBOUR CONDITIONS
                neighbourConditionA = norm(inf(knownObstacles(item).position)) < obj.neighbourDist;     
                if ~neighbourConditionA                                    % If object has any chance of been in the neighbourhood
                    continue
                end
                % OBSTACLE KNOWLEDGE
                r_b = knownObstacles(item).radius;                  
                p_b = knownObstacles(item).position + p_a;        
                v_b = knownObstacles(item).velocity + v_a;     % Convert relative parameters to absolute
                                                                                            
                % COMPUTE THE VECTOR SHARING PROBLEM
                [Voptimal] = obj.defineProblem(desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)
                optimalSet = horzcat(optimalSet,Voptimal);
            end

            % \\\\\\\\\\\\\\ COMPUTE OBSTACLE INTERSECTIONS \\\\\\\\\\\\\\\            
            if ~isempty(optimalSet)
                % GET THE INTESECTION OF THE OPTIMAL REGIONS
                [avoidanceBox] = obj.evaluateAvoidanceIntersections(desiredVelocity,optimalSet);
            else
                % NOTHING TO AVOID, CHOOSE OPTIMAL VELOCITY
                avoidanceBox = desiredVelocity;
            end
            
            % HANDLE SPECIAL CASE - VELOCITY MAGNITUDE IS ZERO
            speed = norm(mid(avoidanceBox));
            headingVector = mid(avoidanceBox)/speed;
            if isnan(headingVector)
                headingVector = [1;0];                                     % Retain previous heading
            end
        end
        % DEFINE THE INTERVAL GEOMETRIC AVOIDANCE PROBLEM %%%%%%%%%%%%%%%%%
        function [U_a] = defineProblem(obj,desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)
            % This function calculates the avoidance vectors based on the
            % principle of vector sharing.
            % INPUTS:
            % desiredVelocity - The true optimal vector
            % p_a - The agent's absolute position
            % v_a - The agent's absolute velocity
            % r_a - The agent's radius
            % p_b - The obstacle's absolute position
            % p_b - The obstacle's absolute velocity
            % p_b - The obstacle's radius
            % visualiseProblem - Plot flag
            % OUTPUTS:
            % U_a - The optimal vector heading, scaled by the desired
            %       velocity.

            % DEFINE THE INITIAL PROBLEM PARAMETERS
            r   = p_b - p_a;                          % The relative position vector = relative seperation
            c_b = v_b - v_a;                          % Define absolute velocity of b
            c_b_unit = obj.intervalUnit(c_b);         % The unit relative velocity of B
           
            % THE 'NEAR-MISS' VECTOR           
            temp = obj.intervalCross([r;0],[c_b_unit;0]);
            r_m  = obj.intervalCross([c_b_unit;0],temp);                       % The near-miss vector
            r_m  = r_m(1:2,1);
                        
            % DEFINE THE TIME TO CLOSEST APPROACH(+ve converging)
            tau = -(dot(r,c_b)/dot(c_b,c_b));
            
            % AVOIDANCE CONDITION
            weThinkWeArePast = mid(tau) < 0;
            isNoChanceOfCollision = sup(tau) < 0;
            isInvalidPrediction = isinf(tau) || isnan(tau);                % [TO CONFIRM]
            if isNoChanceOfCollision || isInvalidPrediction || weThinkWeArePast                                      
               fprintf('avoidance skipped.\n');
               U_a = desiredVelocity;
               return                 
            end
            
            % DEFINE THE REST REGION
            r_safe = r_a + r_b;                                            % The specified minimum safe seperation limit (m)       
            
            
            
            % GET THE NORM OF THE MISS INTERVAL
            rNorm = obj.intervalNorm(r_m) % mag(r_m)?
            

            % DEFINE THE REST REGION      
            r_res  = r_safe - rNorm
            
            % RESOLUTION RATION
            resRatio = r_res/rNorm;
            
            % DEFINE THE VECTOR SHARING TERMS
%             r_vsa = (norm(V_b)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(-r_m)
%             r_vsb = (norm(V_a)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(r_m)
            AShareRatio = obj.intervalNorm(v_b)/(obj.intervalNorm(v_a) + obj.intervalNorm(v_b));  
            r_vsa = AShareRatio*resRatio*-r_m;
                        
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (v_a*tau + r_vsa);
            
            % NORMALISED AND DEMENSION BY DESIRED SPEED
%             U_a = (U_a/obj.intervalNorm(U_a))*norm(desiredVelocity);
        end
        
        % //////////////////// SENSOR FUNCTIONS ///////////////////////////
        % AGENT-SPECIFIC UPDATE FUNCTION (USING SPHERICAL DATA & SENSOR MODEL)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,dt,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = [];
            obstacleSet = [];
            waypointSet = [];   % No observed objects 
            
            % INPUT HANDLING
            if isempty(observedObjects)
                return
            end
            
            % UPDATE THE AGENT'S UNDERSTANDING OF THE ENIVIRONMENT
            for item = 1:length(observedObjects)
                % DEFINE THE SENSED VARIABLES
                % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
                
                % ADDITIONAL MEASUREMENTS
                [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.getSensorMeasurementIntervals(observedObjects(item));
                % CONVERT SPHERICAL MEASUREMENTS TO RELATIVE CARTESIAN
                [measuredPosition] = obj.cartesianFromSpherical(measuredRange,measuredAzimuth,measuredElevation);    % Calculate new relative position
                measuredPosition = measuredPosition(1:2,1); % Reduce 3D inputs to 2D inputs
                            
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                objectRadius = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % STATE ESTIMATION FUNCTION 
                [priorPosition,priorVelocity] = obj.getLastStateFromMemory(observedObjects(item).objectID);
                
                if isempty(priorPosition)
                    % FIRST SIGHT OF OBSTACLE
                    measuredVelocity = NaN(2,1);
                else
                    % CALCULATE THE STATE UPDATE
                    [measuredPosition,measuredVelocity] = obj.linearStateEstimation(dt,priorPosition,priorVelocity,measuredPosition);
                end
                
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(observedObjects(item).priority)
                    objectPriority = observedObjects(item).priority;          % Assign priority to memory
                else
                    objectPriority = 1/observedObjects(item).timeToCollision; % Assign priority based on TTC 
                end
                
                % MEMORY ITEM CONSTRUCTION
                memoryItem = struct('name',observedObjects(item).name,...
                    'objectID',observedObjects(item).objectID,...
                    'type',observedObjects(item).type,...
                    'radius',objectRadius,...
                    'position',measuredPosition,...
                    'velocity',measuredVelocity,...
                    'range',measuredRange,...
                    'azimuth',measuredAzimuth,...
                    'elevation',measuredElevation,...
                    'alpha',measuredAlpha,...
                    'TTC',observedObjects(item).timeToCollision,...
                    'priority',objectPriority); 
                
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.updateAgentKnowledge(memoryItem);
            end
            
            % SORT OBJECT SET BY PRIORITY
            if isintval([obj.memory.priority])
                [~,ind] = sort(mid([obj.memory.priority]),2,'descend');    % Ordered indices of the object IDs
                obj.memory = obj.memory(ind);
            else
                [~,ind] = sort([obj.memory.priority],2,'descend');         % Ordered indices of the object IDs
                obj.memory = obj.memory(ind);
            end
            
            % DECERN OBJECT TYPES
            % obj.memory now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            agentIndices    = [obj.memory.type] == OMAS_objectType.agent;
            waypointIndices = [obj.memory.type] == OMAS_objectType.waypoint;  
            obstacleIndices = [obj.memory.type] == OMAS_objectType.obstacle; % Differenciate between the different object types
            
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet    = obj.memory(agentIndices);
            waypointSet = obj.memory(waypointIndices);
            obstacleSet = obj.memory(obstacleIndices); 
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.waypointUpdater(waypointSet);  % Update the target waypoint and heading vector
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]