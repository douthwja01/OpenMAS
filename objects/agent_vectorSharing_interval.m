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

% %             obj.linearVelocityLimits = [30;30;1];                          % Limits on the agents velocity 
% %             obj.linearAccelerationLimits = [1;1;1];                        % Limits on the agents acceleration
% %             obj.angularVelocityLimits = [pi/6;pi/6;pi/6];
% %             obj.angularAccelerationLimits = [0.1;0.1;0.1]; % rad/s^2
% %             obj.SENSORS
%             % SENSOR UNCERTAINTY
%             obj.positionStandardDeviation = 0.5;          % Accurate to within 0.5m
%             obj.velocityStandardDeviation = 0.1;          % Accurate to within 0.1m/s
%             obj.rangeFinderStandardDeviation = 0.1;      % Accurate to within 0.2m
%             obj.cameraStandardDeviation = 5.208E-5;       % 1 pixel in a 1080p image

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
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
            desiredHeadingVector = [1;0;0];
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
            [newState] = obj.dynamics_simple(dt,desiredVelocity,[0;0;0]);
            obj = obj.updateGlobalProperties_fixedFrame(dt,newState);  
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputNames = {'Vx (m/s)','Vy (m/s)','Vz (m/s)','Pitch (rad)','Roll (rad)','Yaw (rad)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [newState(7:9);newState(4:6)];         % Record the control inputs
            
            % // DISPLAY CONFLICT RESOLUTION
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                obj = obj.getAnimationFrame(overHandle,TIME,'resolutionZone.gif');
                close(overHandle);
            end
        end
    end
    
    methods
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
                [Voptimal_A] = obj.defineProblem(desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)                
                optimalSet = horzcat(optimalSet,Voptimal_A);
            end    
                
            % \\\\\\\\\\\\\\ COMPUTE OBSTACLE INTERSECTIONS \\\\\\\\\\\\\\\            
            if ~isempty(optimalSet)
                % GET THE INTESECTION OF THE OPTIMAL REGIONS
                [avoidanceBox] = obj.evaluateAvoidanceIntersections(desiredVelocity,optimalSet)
%                 [avoidanceBox] = obj.evaluateAvoidanceIntersections_legacy(desiredVelocity,optimalSet)
            else
                % NOTHING TO AVOID, CHOOSE OPTIMAL VELOCITY
                avoidanceBox = desiredVelocity;
            end
                            
            % HANDLE SPECIAL CASE - VELOCITY MAGNITUDE IS ZERO
            speed = mid(obj.intervalNorm(avoidanceBox));
            headingVector = mid(avoidanceBox)/speed;
            if any(isnan(headingVector))
                headingVector = v_a/mid(obj.intervalNorm(v_a));            % Retain previous heading
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
            temp = obj.intervalCross(r,c_b_unit);
            r_m  = obj.intervalCross(c_b_unit,temp);                       % The near-miss vector
            
            % DEFINE THE TIME TO CLOSEST APPROACH(+ve converging)
            tau = -(dot(r,c_b)/dot(c_b,c_b))
            
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
            rNorm = obj.intervalNorm(r_m);% mag(r_m)?
            % DEFINE THE REST REGION      
            r_res  = r_safe - rNorm;
            % RESOLUTION RATION
            resRatio = r_res/rNorm;
            % DEFINE THE VECTOR SHARING TERMS
            % r_vsa = (norm(V_b)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(-r_m)
            % r_vsb = (norm(V_a)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(r_m)
            AShareRatio = obj.intervalNorm(v_b)/(obj.intervalNorm(v_a) + obj.intervalNorm(v_b));  
            r_vsa = AShareRatio*resRatio*-r_m;
                        
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (v_a*tau + r_vsa);
        end
        
        % DEFINE THE INTERVAL GEOMETRIC AVOIDANCE PROBLEM %%%%%%%%%%%%%%%%%
        function [U_a,U_b] = defineProblem_legacy(obj,V_a,S_b,c_b,obstacleSize,visualiseProblem)
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
%             mid(r_m)
            
            tol = 1E-1;
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
            r_safe = obj.VIRTUAL.radius + obstacleSize; % The specified minimum safe seperation limit (m)       
            r_res  = r_safe - rNorm;
            
            % RESOLUTION RATION
            resRatio = r_res/rNorm;
%             resRatio = r_res/norm(r_m)c
            % DEFINE THE VECTOR SHARING TERMS
%             r_vsa = (norm(V_b)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(-r_m)
%             r_vsb = (norm(V_a)/(norm(V_a) + norm(V_b)))*(r_res/norm(r_m))*(r_m)
            AShareRatio = obj.intervalNorm(V_b)/(obj.intervalNorm(V_a) + obj.intervalNorm(V_b));
            BShareRatio = obj.intervalNorm(V_a)/(obj.intervalNorm(V_a) + obj.intervalNorm(V_b));
                        
            r_vsa = AShareRatio*resRatio*-r_m;
            r_vsb = BShareRatio*resRatio*r_m;
                        
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (V_a*tau + r_vsa);
            unit_Ua = obj.intervalUnit(U_a);
            U_b = (V_b*tau + r_vsb);
            unit_Ub = obj.intervalUnit(U_b);
            
            % DEFINE THE OPTIMAL ACCELERATION VECTOR
            % Based on the hamiltonian, H, the vector that minimises the
            % cost expression J, is an ahat that is parallel to -r_m:
            % J = -0.5*norm(r_m)^2;
            % H = dot(-r_m,c_b) + (a_mod*tau)*dot(-r_m,unit_acc);
            
            % We can therefore write the optimal acceleration direction:
            % unit_A = -obj.intervalUnit(r_m);
        end
    end
    % //////////////////////// SENSOR FUNCTIONS ///////////////////////////     
    methods
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
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL 
                objectRadius = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % STATE ESTIMATION FUNCTION 
                [priorPosition,priorVelocity] = obj.getLastStateFromMemory(observedObjects(item).objectID);
                
                if isempty(priorPosition)
                    % FIRST SIGHT OF OBSTACLE
                    measuredVelocity = NaN(3,1);
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
        % SENSOR MODEL - LOCAL IMU & GPS INTERVALS (2D & 3D)
        function [position,velocity,radius] = getAgentMeasurementIntervals(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).

            if obj.VIRTUAL.is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices = 1:2;
                velocityIndices = 4:5;
            end
            % GET THE STATE UNCERTAINTIES
            positionAUncertainty = obj.SENSORS.sigma_position*randn(numel(positionIndicies),1);  % Zero estimate of position
            velocityAUncertainty = obj.SENSORS.sigma_velocity*randn(numel(velocityIndicies),1);  % Zero estimate of velocity
            % GET THE LOCAL ABSOLUTE VARIABLES
            position = obj.localState(positionIndices,1) + midrad(positionAUncertainty,3*obj.SENSORS.sigma_position);
            velocity = obj.localState(velocityIndices,1) + midrad(velocityAUncertainty,3*obj.SENSORS.sigma_velocity); 
            % RADIUS IS ASSUMED KNOWN
            radius = obj.VIRTUAL.radius;
        end
        % SENSOR MODEL - CAMERA & RANGE FINDER INTERVALS (2D & 3D)
        function [rangeBox,azimuthBox,elevationBox,alphaBox] = getSensorMeasurementIntervals(obj,obstacleData) 
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
            [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.getSensorMeasurements(obstacleData);
            % DEFINE THE MEASUREMENT INTERVALS
            rangeBox     = midrad(measuredRange,obj.confidenceAssumption*obj.SENSORS.sigma_rangeFinder);
            azimuthBox   = midrad(measuredAzimuth,obj.confidenceAssumption*obj.SENSORS.sigma_camera);
            elevationBox = midrad(measuredElevation,obj.confidenceAssumption*obj.SENSORS.sigma_camera);
            alphaBox     = midrad(measuredAlpha,obj.confidenceAssumption*obj.SENSORS.sigma_camera);
        end
    end
%   STATIC METHODS    
    methods (Static) % Don't require an input object to work        
        % CALCULATE THE GLOBAL OPTIMAL INTERVAL
        function [optimalBox] = evaluateAvoidanceIntersections(desiredVelocity,avoidanceVelocities)
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
               optimalBox = desiredVelocity;
               return
            end
                        
            % MOVE THROUGH THE VELOCITY SET
            viableIntersections = avoidanceVelocities(:,1);                          % Initialise the global intersection with the highest priority 
            
            for i = 1:(size(avoidanceVelocities,2))% + 1)
                % GET THE CANDIDATE VELOCITY
                if i == (size(avoidanceVelocities,2) + 1)
                    candidateInterval = desiredVelocity;
                else
                    candidateInterval = avoidanceVelocities(:,i);
                end
                
                if isnan(candidateInterval)
                    warning('Velocity candidate is NaN');
                end
                
                % THE EMPTY INTERSECTION
                [failureFlags,temp] = emptyintersect(viableIntersections,candidateInterval);
                
                if sum(failureFlags) ~= 0
                    % Some elements do not successfully intersect \|~|/  
                else
                    % The regions intersected successfully
                    viableIntersections = temp;
                end
            end
            
            % DEFINE THE INTERSECTION AS THE LOCALLY OPTIMAL AVOIDANCE SET
            optimalBox = viableIntersections;
        end
        % CALCULATE THE GLOBAL OPTIMAL INTERVAL
        function [optimalHeading] = evaluateAvoidanceIntersections_legacy(waypointHeading,avoidanceVelocities)
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
        end
    end
    % INTERVAL METHODS & MATH TOOLS
    methods (Static)
        % CONVERT CARTESIAN TO SPHERICAL
        function [range,azimuth,elevation] = sphericalFromCartesian(positionVector)
            range = norm(positionVector);                                                           % The range
            azimuth = atan(positionVector(2)/positionVector(1));                                    % The angle made in the azimuth (bearing)
            elevation = atan(positionVector(3)/sqrt(positionVector(1).^2 + positionVector(2).^2));  % The elevation (vertical bearing)
        end
        % INTERVAL CROSS PRODUCT
        function [crossInt] = intervalCross(intA,intB)
            crossInt = [intA(2)*intB(3)- intA(3).*intB(2); 
                        intA(3)*intB(1)- intA(1).*intB(3); 
                        intA(1)*intB(2)- intA(2).*intB(1)]; 
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