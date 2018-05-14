%% THE 2D AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the 'agent' class, this object class is intended provide
% support for 2D simulation. The functions provided here 

% Author: James A. Douthwaite

classdef agent_2D < agent
%%% AGENT(2D) BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % PROPERTIES UNIQUE TO THE RECIPROCAL VO METHOD
    end
%%   CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_2D(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin); % Call the super class
            % THE 2D STATE VECTOR
            obj.localState = zeros(6,1);
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % ///////////////// AGENT MAIN CYCLE //////////////////////////////
        function obj = processTimeCycle(obj,TIME,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % DEFAULT BEHAVIOUR 
            desiredSpeed = 2;
            desiredHeadingVector = [1;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
                desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
            end
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            headingRate = 0.00;      % Constant +ve yaw rate
            speed = 1;              % Constant speed (m/s)
            
            if TIME.currentTime >= 0 && TIME.currentTime < 5
                headingRate = 0.5;
                fprintf('yawing...\n');
            end
            
            if TIME.currentTime >= 5 && TIME.currentTime < 10
                headingRate = -0.5;
                fprintf('pitching...\n');
            end
            
            if TIME.currentTime >= 10 && TIME.currentTime < 15
                headingRate = 0.5;
                fprintf('rolling...\n');
            end   
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_singleIntegrator(dt,[speed;0],headingRate);
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
    end
    
    % ///////////////////// (2D) SENSORY FUNCTIONS ////////////////////////
    methods
        % ALL IN ONE AGENT-INFORMATION UPDATE (FOR 2D STATES)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment in 2D.
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
            
            % ONLY TAKE THE 2D STATES
            linearIndicies = 1:2;
            
            for item = 1:length(observedObjects)
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(observedObjects(item).priority)
                    objectPriority = observedObjects(item).priority;           % Assign priority to memory
                else
                    objectPriority = 1/observedObjects(item).timeToCollision; 
                end
                
                % MEMORY ITEM CONSTRUCTIONS
                memoryItem = struct('name',observedObjects(item).name,...
                                'objectID',observedObjects(item).objectID,...
                                    'type',observedObjects(item).type,...
                                  'radius',observedObjects(item).radius,...
                                'position',observedObjects(item).position(linearIndicies,1),...
                                'velocity',observedObjects(item).velocity(linearIndicies,1),...
                                     'TTC',observedObjects(item).timeToCollision,...
                                'priority',objectPriority); 
                
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.updateAgentKnowledge(memoryItem);
            end
            
            % SORT OBJECT SET BY PRIORITY
            [~,ind] = sort([obj.memory.priority],2,'descend');             % Ordered indices of the object IDs
            obj.memory = obj.memory(ind);
            
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
            [obj,~] = obj.waypointUpdater(waypointSet);                    % Update the target waypoint and heading vector
        end
        % AGENT-SPECIFIC UPDATE FUNCTION (USING SPHERICAL DATA & SENSOR MODEL)
        function [obj,obstacleSet,agentSet,waypointSet] = getSensorUpdate(obj,dt,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observationSet - The full observed object structure
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
                [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.getSensorMeasurments(observedObjects(item));
                [measuredPosition] = obj.cartesianFromSpherical(measuredRange,measuredAzimuth,measuredElevation);    % Calculate new relative position
                % MAP TO MEASURED POSITION
                measuredPosition = measuredPosition(1:2,1);
                
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                measuredRadius = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % STATE ESTIMATION FUNCTION 
                [priorPosition,priorVelocity] = obj.getLastStateFromMemory(observedObjects(item).objectID);       % Get prior obstacle knowledge
                
                if isempty(priorPosition)
                    % FIRST SIGHT OF OBSTACLE
                    measuredVelocity = NaN(2,1);
                else
                    % CALCULATE THE STATE UPDATE
                    [measuredPosition,measuredVelocity] = obj.linearStateEstimation(dt,...
                    priorPosition,priorVelocity,measuredPosition);
                end
                                
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(observedObjects(item).priority)
                    objectPriority = observedObjects(item).priority;       % Assign priority to memory
                else
                    objectPriority = 1/observedObjects(item).timeToCollision; 
                end
                
                % MEMORY ITEM CONSTRUCTION
                memoryItem = struct('name',observedObjects(item).name,...
                    'objectID',observedObjects(item).objectID,...
                    'type',observedObjects(item).type,...
                    'radius',measuredRadius,...
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
            [~,ind] = sort([obj.memory.priority],2,'descend');             % Ordered indices of the object IDs
            obj.memory = obj.memory(ind);

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
        
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [position,velocity,radius] = getAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            positionIndices = 1:2;
            velocityIndices = 4:5;
            % GET THE LOCAL ABSOLUTE VARIABLES
            position =  obj.localState(positionIndices,1) + obj.SENSORS.positionSigma*rand(2,1); % Get the absolute position measurement
            velocity =  obj.localState(velocityIndices,1) + obj.SENSORS.velocitySigma*rand(2,1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            radius = obj.VIRTUAL.radius;
        end
    end
    
    % //////////////////////// DYNAMICS & CONTROL //////////////////////////
    methods   
        % IMPOSE 2D KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = kinematicContraints(obj,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % HEADING IS RELATIVE
            currentVelocity = obj.localState(4:5,1);
            currentHeading = zeros(size(desiredHeading));
            
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - currentHeading;
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = obj.boundValue(Hdot,-obj.DYNAMICS.angularVelocityLimits,obj.DYNAMICS.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = currentHeading + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed;                          % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = obj.boundValue(acceleration,-obj.DYNAMICS.linearAccelerationLimits(1),obj.DYNAMICS.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = obj.boundValue(achievedSpeed,-obj.DYNAMICS.linearVelocityLimits(1),obj.DYNAMICS.linearVelocityLimits(1)); 
        end
        % INITIALISE 2D DYNAMIC PARAMETERS
        function [obj] = getDynamicParameters(obj)
            % This function calls the dynamics structure that gives the
            % agent dynamic/kinematic limits and parameters.
            % BUILD THE DYNAMICS SUB-STRUCTURE
            obj.DYNAMICS = struct('linearVelocityLimits',[inf;inf],...     % Limits on the agents linear velocity 
                              'linearAccelerationLimits',[inf;inf],...     % Limits on the agents linear acceleration
                                 'angularVelocityLimits',[inf;inf],...     % Limits on the agents angular velocity
                             'angularAccelerationLimits',[inf;inf]);       % Limits on the agents angular acceleration
        end
    end
    
    % //////////////////// SIMULATION & CORE INTERFACES ///////////////////
    methods
        % THE DEFAULT 2D STATE INITIALISER
        function [obj] = initialise_localState(obj,localENUVelocity,localENUrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y phi dx dy dphi]
            
            % Build the initial state vec
            obj.localState = zeros(6,1);
            obj.localState(3) = localENUrotations(3);     % The initial heading angle
            obj.localState(4:5) = localENUVelocity(1:2);  % The initial 2D velocities
        end
        % UPDATE THE GLOBAL PROPERTIES ( LOCAL AXIS REMAINS STATIC )
        function [obj] = updateGlobalProperties_fixedFrame(obj,dt,localENUState)
            % This function computes the new global parameters for the given
            % agent based on its new state.
            
            % CONSTANT PARAMETERS
            linearIndices = 4:5;
            
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                obj.VIRTUAL.idleStatus = 1;
                localENUState(4:6) = zeros(3,1); % Freeze the agent
            end
            
            % //////////////// DEFINE PREVIOUS PARAMETERS /////////////////
            quaternion_k        = obj.VIRTUAL.quaternion;
            velocity_ENU_k_plus = localENUState(linearIndices,1);
            
            % ////////////////// [WORKING: NON-ROTATED] ///////////////////
            quaternion_k_plus = quaternion_k;
            
            %NEW ROTATION MATRIX
            [R_update] = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            % THE VELOCITY VECTORS AT K & K+1
            localVelocityUpdate = [velocity_ENU_k_plus;0];                 % 3D although locally 2D
            
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = R_update*localVelocityUpdate;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*globalVelocity_k_plus;
            
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.rotationMatrix = R_update;
            obj.localState = localENUState;
        end
        % UPDATE THE GLOBAL REPRESENTATION ( FOR 2D AGENTS )    [CONFIRMED]
        function [obj] = updateGlobalProperties(obj,dt,localNEDState)
            % This function updates the agents global properties from
            % X_global(@t=k] to X_global(@t=k+1) for a 2D agent defined by
            % state vector: [x y lambda dx dy dlambda]
            
            % CONSTANT PARAMETERS
            angularIndices = 3;
            linearIndices = 4:5;
            
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                obj.VIRTUAL.idleStatus = 1;
                localNEDState(4:6) = zeros(3,1); % Freeze the agent
            end
            
            % CALCULATE THE AXIS RATES FROM THE NEW STATE
            yawRate = (localNEDState(angularIndices) - obj.localState(angularIndices))/dt;
            
            % //////// GET THE UPDATED QUATERNION /////////////////////////
            globalRates = zeros(3,1);
            globalRates(1) = -yawRate; % rates about the Z axis respectively
            
            [quaternion_k_plus] = OMAS_axisTools.updateQuaternion(obj.VIRTUAL.quaternion,globalRates,dt);
            R_new = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = R_new*[localNEDState(linearIndices);0];
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*globalVelocity_k_plus;
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.priorState = obj.localState;
            obj.localState = localNEDState;                                % Reassign the obj.localstate
        end
    end
end
