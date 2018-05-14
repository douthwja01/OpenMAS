%% THE AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic agent and import this variables 
% into the simulation space for the purpose of multi-vehicle control simulation.
% The agent object is a child of the objectDefintion; the prinicple
% distinctions being:
% sensorRange      - The agent is assumed capable of observing its
%                    surroundings.
% controlFrequency - The frequency at which the control cycle is computed.

% Author: James A. Douthwaite

classdef agent < objectDefinition
%%% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % OBSTACLE OBJECT MATRIX (record of sightings)
        memory;                         % Record of last known [objectID;positions;velocities]
        % WAYPOINTS
        targetWaypoint;                 % The current waypoint target
        achievedWaypoints;              % The agents list of locally achieved waypoints
        % DYNAMIC PARAMETERS
        DYNAMICS;
        % SENSOR MEASUREMENT PARAMETERS
        SENSORS;
        % AGENT-SIDE OUTPUT DATA
        DATA;                           % The output container for agent-side data.
    end
%%   CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@objectDefinition(varargin); % Call the super class
            
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            
            % THE 2D STATE VECTOR
            obj.localState = zeros(12,1);
                                    
            % VIRTUAL DEFINITION
            obj.VIRTUAL.type = OMAS_objectType.agent;
            obj.VIRTUAL.symbol = 'diamond';
            obj.VIRTUAL.detectionRange = inf;  
            obj.VIRTUAL.idleStatus = logical(false);
            % Assume the agent has perfect environmental knowledge (m)

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % ///////////////////// AGENT MAIN CYCLE //////////////////////////
        function obj = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
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
            desiredHeadingVector = [1;0;0];
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
            headingRate = 0.1;      % Constant +ve yaw rate
            speed = 1;              % Constant speed (m/s)
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_singleIntegrator(dt,[speed;0;0],[0;0;headingRate]);
            
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
    end
    % ////////////////////// (3D) SENSORY FUNCTIONS ///////////////////////
    methods
        % ALL IN ONE AGENT-INFORMATION UPDATE (FROM ENVIRONMENT)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,observedObjects)
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
                                'position',observedObjects(item).position,...
                                'velocity',observedObjects(item).velocity,...
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
                
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                measuredRadius = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % STATE ESTIMATION FUNCTION 
                [priorPosition,priorVelocity] = obj.getLastStateFromMemory(observedObjects(item).objectID);       % Get prior obstacle knowledge
                
                if isempty(priorPosition)
                    % FIRST SIGHT OF OBSTACLE
                    measuredVelocity = NaN(3,1);
                else
                    % CALCULATE THE STATE UPDATE
                    [measuredPosition,measuredVelocity] = obj.linearStateEstimation(dt,priorPosition,priorVelocity,measuredPosition);
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
        
        % SENSOR MODEL - CAMERA & RANGE FINDER
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
            range     = obstacleData.range + obj.SENSORS.rangeFinderSigma*randn(1);     
            azimuth   = obstacleData.azimuthAngle + obj.SENSORS.cameraSigma*randn(1);
            elevation = obstacleData.inclinationAngle + obj.SENSORS.cameraSigma*randn(1);
            % OBSERVED RADIUS
            angularWidth = obstacleData.angularWidth + obj.SENSORS.cameraSigma*randn(1);      % Uncertainty in the angular measurement
        end
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [position,velocity,radius] = getAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            positionIndices = 1:3;
            velocityIndices = 7:9;
            % GET THE LOCAL ABSOLUTE VARIABLES
            position =  obj.localState(positionIndices,1) + obj.SENSORS.positionSigma*rand(3,1); % Get the absolute position measurement
            velocity =  obj.localState(velocityIndices,1) + obj.SENSORS.velocitySigma*rand(3,1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            radius = obj.VIRTUAL.radius;
        end
        % SENSOR MODEL - PERFECT SENSING
        function [obj] = getSensorParameters(obj)
            % This function is designed to populate the SENSOR field with
            % representative sensor uncertainty.
            % BUILD THE SENSOR FACTOR
            obj.SENSORS = struct('sensorHorizon',inf,...    % Assume the agent has perfect environmental knowledge (m)
                                 'positionSigma',0.0,...    % Accurate to within 0.5m
                                 'velocitySigma',0.0,...    % Accurate to within 0.1m/s
                              'rangeFinderSigma',0.0,...    % Accurate to within 0.1m
                                   'cameraSigma',0.0,...    % One pixel in a 1080p image
                               'sampleFrequency',inf);      % Object has perfect precision
        end
    end
    
    % //////////////////////// DATA I/O FUNCTIONS /////////////////////////
    methods
        % AGENT DATA STORAGE
        function [obj] = writeAgentData(obj,TIME,loopIndicator,loopDuration)
            % This function writes a regular DATA structure to the agent
            % .DATA structure to allow it to be easily interpretted.
                        
            % THE AGENT COMPUTATION TIMES
            obj.DATA.algorithm_indicator(TIME.currentStep) = loopIndicator;
            obj.DATA.algorithm_dt(TIME.currentStep) = loopDuration;
            obj.DATA.algorithm_steps(TIME.currentStep) = TIME.currentStep;
            obj.DATA.algorithm_time(TIME.currentStep) = TIME.currentTime;
        end 
    end
    
    % //////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods    
        % SPEED AND HEADING (4D) PID CONTROLLER
        function [d_speed,d_heading,obj] = trajectoryController(obj,targetVelocity)
            % This function is desiged to generate feedack from a desired
            % local velocity vector. 
            % INPUTS:
            % targetHeading - The unit heading vector
            % targetSpeed   - The target speed in that direction
            % OUTPUTS:
            % heading_fb    - The feedback on the agent heading
            % speed_fb      - The feedback on the agent speed
            % obj           - The error-updated agent object
            
            % VELOCITY PARAMETERS
            targetSpeed = norm(targetVelocity);
            targetHeading = targetVelocity/targetSpeed;
            
            % DEFINE THE STATE ERROR TERMS
            % RELATIVE HEADING
            [e_theta,e_lambda] = obj.getControlInputs([1;0;0],targetHeading);  % Get the equivalent angular heading and elevation
            % RELATIVE SPEED
            e_speed = targetSpeed - norm(obj.localState(7:9));                 % The speed error 
            controlError = [e_speed;0;e_theta;e_lambda];                       % -ve in the NED control frame of reference
            
            % TUNING PARAMETERS 
            Kp_linear = 0.8;
            Kd_linear = 0;
            Kp_angular = 1;
            Kd_angular = 0;
            
            % CALCULATE THE CONTROL FEEDBACK
            control_fb = diag([Kp_linear Kp_angular Kp_angular Kp_angular])*controlError + ...
                         diag([Kd_linear Kd_angular Kd_angular Kd_angular])*(controlError - obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = controlError;
            
            % PARSE THE CONTROL FEEDBACK FOR EXTERNAL USE
            d_speed    = control_fb(1);             % Absolute speed input
            d_heading  = control_fb(2:4);
            % CONVERT FLU (AVOIDANCE TO NED (CONTROL)
            d_heading = diag([1 1 -1])*d_heading;
        end       
        % SIMPLE CONTROLLER
        function [d_speed,d_heading] = simpleController(obj,targetVelocity)
            % This function preforms the most basic of update functions. We
            % simply translate the desired velocity (in the body axis) into
            % the new NED state vector at time k.
            
            % ASSUMPTIONS:
            % - The agent is orientated with the body x axis and velocity 
            %   vector colinear (and otherwise S&L).
            % - The velocity vector is provided in the body frame of
            %   reference.
            
            % OUTPUTS
            % d_heading - Feedback on the vector angles.
            % d_speed   - Feedback on the speed. 
            

            % INPUT HANDLING
            inputDim = numel(targetVelocity);
            referenceAxis = zeros(size(targetVelocity));
            referenceAxis(1) = 1;
            
            % VELOCITY PARAMETERS
            targetSpeed = norm(targetVelocity);
            targetHeading = targetVelocity/targetSpeed;
            
            % DECOMPOSE THE HEADING FEEDBACK FROM THE ANGULAR STATES
            if inputDim == 2
                d_psi = sign(targetHeading(2))*acos(dot(referenceAxis,targetHeading)/norm(referenceAxis)); % Only yaw is altered
                d_heading = -d_psi;            
                % COMPUTE FORWARD SPEED FEEDBACK
                d_speed = targetSpeed - norm(obj.localState(4:5,1));                                        % The heading change is yaw only
            else
                [d_psi,d_theta] = obj.getVectorHeadingAngles(referenceAxis,targetHeading); 
                d_heading = [0;d_theta;-d_psi];                            % Roll is undefined from the velocity alone
                % COMPUTE FORWARD SPEED FEEDBACK
                d_speed = targetSpeed - norm(obj.localState(7:9,1));
            end
        end
        % IMPOSE KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = kinematicContraints(obj,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            currentVelocity = obj.localState(7:9);
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - zeros(3,1);
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = obj.boundValue(Hdot,-obj.angularVelocityLimits,obj.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = zeros(3,1) + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed; % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = obj.boundValue(acceleration,-obj.linearAccelerationLimits(1),obj.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = obj.boundValue(achievedSpeed,-obj.linearVelocityLimits(1),obj.linearVelocityLimits(1)); 
        end
        % INITIALISE DYNAMICS
        function [obj] = getDynamicParameters(obj)
            % This function calls the dynamics structure that gives the
            % agent dynamic/kinematic limits and parameters.
            % BUILD THE DYNAMICS SUB-STRUCTURE
            obj.DYNAMICS = struct('linearVelocityLimits',[inf;inf;inf],... % Limits on the agents linear velocity 
                              'linearAccelerationLimits',[inf;inf;inf],... % Limits on the agents linear acceleration
                                 'angularVelocityLimits',[inf;inf;inf],... % Limits on the agents angular velocity
                             'angularAccelerationLimits',[inf;inf;inf]);   % Limits on the agents angular acceleration
        end
    end
    
    % ////////////////////// GENERAL STATIC METHODS ///////////////////////
    methods (Static)
        % DEFINE THE PITCH AND YAW TO MEET TARGET HEADING
        function [lambda,theta] = getVectorHeadingAngles(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS:
            % V - The current unity velocity vector
            % U - The unit correction vector
            % OUTPUTS:
            % lambda - The azimuth angle (2D)
            % theta  - The pitch angle   (3D)
                        
            % GET THE LINE OF SIGHT ANGLE             
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];                                            % Reject the vertical elements
            % DEFINE LOS ANGLE            
            rotationAxis = cross(Vh,Uh);
            lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            % GET THE ELEVATION ANGLE 
            theta = atan2(U(3),norm(Uh));
            
            % CHECK FOR PARALLISM
%             if sign(dot(V,U)) == -1
% %                 display('opposing vector');
%                 lambda = pi - lambda;
%                 theta = pi - theta;
%             end
        end
        % CALCULATE THE NEW STATE ESTIMATE
        function [newPosition,newVelocity] = linearStateEstimation(dt,position0,velocity0,newPosition)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
                      
           % GET THE POSITION
           dX = (newPosition - position0);
           newVelocity = dX/dt;  % Defines the average velocity
           
           % NO PREVIOUS VELOCITY RECORDED
           if any(isnan(velocity0))
              newPosition = newPosition;
              newVelocity = newVelocity;
              return
           end          
        end
        % GET ESTIMATE OF RELATIVE POSITION & VELOCITY
        function [cartesianPosition] = cartesianFromSpherical(range,azimuth,elevation)
           % This function calculates the relative position from a spherical 
           % coordinate system that assumes angles are measured from the 
           % agents heading vector Xref.
           % INPUTS:
           % range     - The objects radial seperation (m)
           % azimuth   - The objects relative azimuth angle (rad)
           % elevation - The objects relative elevation (rad)
           % OUTPUTS:
           % cartesianPosition  - The new 3D position in local coordinates (m)

           % DEFINE THE POSITION AS A VECTOR INTERVAL
           cartesianPosition = [cos(azimuth)*cos(elevation);...
                                sin(azimuth)*cos(elevation);...
                                             sin(elevation)]*range;
        end
        % CONVERT CARTESIAN TO SPHERICAL
        function [range,azimuth,elevation] = sphericalFromCartesian(positionVector)
            range = norm(positionVector);                                            % The range
            azimuth = atan2(positionVector(2),positionVector(1));                    % The angle made in the azimuth (bearing)
            elevation = atan2(positionVector(3),sqrt(positionVector(1).^2 + positionVector(2).^2));  % The elevation (vertical bearing)
        end
        % VALIDATE THE OBSTACLE
        function [tau] = validateCollision(S_b,c_b)
            % CONFIRM WE KNOW THE HEADING OF THE OBSTACLE
            if any(isnan(c_b))
                tau = -inf;
                return
            end                 
            % Define the time to closest approach (+ve converging)
            tau = -(dot(S_b,c_b)/dot(c_b,c_b));
        end
    end  
    
    % //////////////////// SIMULATION & CORE INTERFACES ///////////////////
    methods
        % SELECT WAYPOINT
        function [obj,waypointVector] = waypointUpdater(obj,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if ~exist('waypointSet','var') || isempty(waypointSet)
                obj.targetWaypoint = [];       % No waypoints are available
                waypointVector = [];           % Default heading
                return
            end
            
            % WAYPOINT SET IS POPULATED -> GET THE PRIORITY VECTOR
            priorityVector = [waypointSet.priority];                       % The vector of priorities
            waypointIndex = linspace(1,length(priorityVector),length(priorityVector)); % Memory of their position in the waypoint set
            waypointIDset  = [waypointSet.objectID];                       % The vector of IDs 
            waypointMatrix = [waypointIndex;waypointIDset;priorityVector];
            
            % WAYPOINT MATRIX IS OF THE FORM
            % [ 1 2 3 ] - Positions in the waypoint set.
            % [ 2 7 4 ] - Waypoint IDs.
            % [ 3 5 9 ] - Waypoint priorities.
                        
            % TARGET WAYPOINT IS EMPTY -> POPULATE WITH NEW WAYPOINTS
            if isempty(obj.targetWaypoint)
                % NO WAYPOINTS HAVE BEEN SELECTED
                if isempty(obj.achievedWaypoints)
                    [~,maxPriorityIndex] = max(waypointMatrix(3,:));       % Index of maximum priority value in waypoint matrix
                    waypointSetIndex = waypointMatrix(1,maxPriorityIndex); % Gets the associated waypoint set index of the highest priority
                    obj.targetWaypoint = waypointSet(waypointSetIndex);    % Update with the highest priority waypoint
                else
                    % SELECT HIGHEST PRIORITY, BUT NOT ACHIEVED
                    invalidWaypointIDs = obj.achievedWaypoints;            % Get the achieved waypoint IDs
                    availableWaypointIDs = waypointMatrix(2,:);            % Get those visable IDs
                    validWaypoints = availableWaypointIDs ~= invalidWaypointIDs;   % Get those visible IDs that are valid
                    % IF NO FURTHER VALID IDs
                    if ~any(validWaypoints)
                        obj.targetWaypoint = [];
                        waypointVector = [];
                        return
                    else
                        % SELECT REMAINING WAYPOINTS
                        validWaypoints = waypointMatrix(:,validWaypoints); % Get waypoinSet location, IDs and priority of valid waypoints
                        % SELECT WAYPOINT OF NEXT HIGHEST PRIORITY                        
                        [~,maxValidIndex ] = max(validWaypoints(3,:));       % Get the index of max priority waypoint
                        waypointSetIndex = validWaypoints(1,maxValidIndex);  % Get the location of the target waypoint in the waypoint set
                        obj.targetWaypoint = waypointSet(waypointSetIndex);  % Select waypoint object
                    end
                end
            end
            
            % EVALUATE CURRENT TARGET WAYPOINT ////////////////////////////          
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            waypointAchievedCondition = norm(obj.targetWaypoint.position) - (obj.targetWaypoint.radius + obj.VIRTUAL.radius);
            % IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
            if waypointAchievedCondition <= 0 && ~any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID))
                obj.achievedWaypoints = horzcat(obj.achievedWaypoints,obj.targetWaypoint.objectID); % Add achieved objectID to set
            end
            
            % CHECK IF CURRENT TARGET IS VALID ////////////////////////////
            % ASSESS THE ACHIEVED WAYPOINT SET AGAINST THE CURRENT TARGET
            invalidTargetCondition = any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID)); % Is the current target ID in the achieved set
            if invalidTargetCondition
               % TARGET WAYPOINT HAS BEEN ACHIEVED, REALLOCATE ONE
               selectionVector = priorityVector < obj.targetWaypoint.priority;   % Only those with lower priority
               if ~any(selectionVector)
                   obj.targetWaypoint = [];
                   waypointVector = [];
                   return
               else
                   reducedMatrix = waypointMatrix(:,selectionVector);            % Build vector of priorities less than 
                   [~,reducedIndex] = max(reducedMatrix(3,:));                   % Finds the max priority index where less than current
                   waypointIndex = reducedMatrix(1,reducedIndex);                % Get the waypoint index from the reduced priority matrix               
                   obj.targetWaypoint = waypointSet(waypointIndex);              % Select the next highest priority waypoint
               end
            else          
                currentWaypointID = obj.targetWaypoint.objectID;                 % Get the ID of the current target waypoint
                selector = (waypointMatrix(2,:) == currentWaypointID);           % Find where the ID appears in the waypoint matrix
                waypointIndex = waypointMatrix(1,selector);                      % Get the waypoint Set index from the waypoint matrix
                obj.targetWaypoint = waypointSet(waypointIndex);                 % Update with the highest priority waypoint  
            end           
            % THE CORRECT TARGET WAYPOINT IS NOW DEFINED, GENERATE HEADING                      
            waypointVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);  
        end
        % UPDATE OBSTACLE KNOWLEDGE
        function [obj] = updateAgentKnowledge(obj,newEntry)
            % Updates the agents knowledge based on the new information
            % sensed by the agent. 
            % INPUTS:
            % objectID   - The objects ID
            % objectType - The objects sim type (obstacle or waypoint)
            % radius     - The objects apparent size
            % newState   - The objects estimated state
            % priority   - The objects priority
            
            % INPUT HANDLING
            if ~isfield(newEntry,'priority')
                newEntry.priority = 0;
            end

            % CHECK CURRENT KNOWLEDGE STATUS
            if isempty(obj.memory)
                obj.memory = vertcat(obj.memory,newEntry); % No knowledge at all
                return
            end
            
            % GET THE LOGICAL INDICES OF ID OCCURANCES
            logicalIDIndex = [obj.memory.objectID] == newEntry.objectID;   % Appearance of object ID in memory
            IDoccurances = sum(logicalIDIndex);
            
            % CHECK FOR THE APPEARANCE OF THE OBJECT IN MEMORY
            if IDoccurances > 1
                error('MEMORY DUPLICATE DETECTED.');
            elseif IDoccurances == 0
                % APPEND ENTRY FOR NEW OBJECT
                obj.memory = vertcat(obj.memory,newEntry); % Append to knowledge to other records
            else
                % UPDATE THE MEMORY STRUCTURE WITH THE NEW ENTRY
                obj.memory(logicalIDIndex) = newEntry;                     % Update existing knowledge of  that agent
            end
        end
        % GET OBJECT MEMORY ITEM FOR OBJECT-ID 
        function [lastMemoryItem] = getLastEntryFromMemory(obj,objectID)
            % This function is designed to return the complete memory item
            % for a given object ID
            % INPUT:
            % obj.memory - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % lastMemoryItem - The agents last memory record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.memory) || ~any([obj.memory.objectID] == objectID)
                lastMemoryItem = [];
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.memory.objectID] == objectID;
            lastMemoryItem = obj.memory(IDindex); 
        end
        % GET STATE FOR OBSTACLE ID 
        function [position,velocity] = getLastStateFromMemory(obj,objectID)
            % This function is designed to return the last known state
            % vector for given obstacle object ID.
            % INPUT:
            % obj.memory - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % lastState - The agents last state record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.memory) || ~any([obj.memory.objectID] == objectID)
                position = [];
                velocity = [];
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.memory.objectID] == objectID;
            position = obj.memory(IDindex).position; 
            velocity = obj.memory(IDindex).velocity; 
        end
        
        % INTIALISE LOCAL 12DOF NED STATE VECTOR 
        function [obj] = initialise_localState(obj,localENUVelocity,localENUrotations)
            % This function is called in order to build a local state
            % vector intended for aerospace control simulation (NED).
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
            
            % BUILD THE DEFAULT ENU STATE VECTOR
%             localFLUState = zeros(12,1);
%             localFLUState(4:6) = localENUrotations;                        % Get the initial local euler heading
%             localFLUState(7:9) = localENUVelocity;                         % Get the initial local velocity vector
                        
            % CONVERT THE FLU (LOCAL ENU) STATES TO THE FRD (LOCAL NED)
            localFRDState = zeros(12,1);
            localFRDState(4:6) = diag([1 1 -1])*localENUrotations;
            
            % ROTATE THE ENU VELOCITY VECTOR TO AS IT APPEARS IN THE NED
            % FRAME
            [R_BG,~] = OMAS_axisTools.eulerToRotationMatrix([(localENUrotations(1)+pi);0;0]);
            localFRDState(7:9) = R_BG*localENUVelocity;
            
            % ASSIGN THE LOCAL FRD STATE
            obj.VIRTUAL.priorState = localFRDState;
            obj.localState = localFRDState;
        end
        % THE GLOBAL PROPERTY FOR GENERAL OBJECT CLASS DERIVATIVES
        function [obj] = updateGlobalProperties(obj,dt,localNEDState)
            % This function is used to calculate the new global (.VIRTUAL)
            % parameters for the current object.

            % updates the agents global properties from
            % X_global(@t=k] to X_global(@t=k+1)
            
            % globalPosition - The global ENU position
            % globalVelocity - The global ENU velocity
            % quaternion - Rotation between the local ENU and global ENU
            
            % CONSTANT PARAMETERS
            velocityIndices = 7:9;
            omegaIndices = 10:12;

            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                obj.VIRTUAL.idleStatus = 1;
                localNEDState(7:12) = zeros(6,1); % Freeze the agent
            end   
            
            % DEFINE UPDATE PARAMETERS
%             ENUrates = localENUState(10:12); % rates about the X Y Z respectively
%             [quaternion_k_plus] = OMAS_axisTools.updateQuaternion(obj.VIRTUAL.quaternion,ENUrates,dt)
%             R_ENU = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);

            % [ TO BE CONFIRMED ] ENU ROTATIONS ARE CURRENTLY MAPPED TO THE
            % 'body rates' - globalRates 
            NEDrates = localNEDState(omegaIndices); % rates about the X Y Z respectively
            globalRates = zeros(3,1);
            globalRates(1) = -NEDrates(3);
            globalRates(2) = NEDrates(2);
            globalRates(3) = NEDrates(1);
            [quaternion_k_plus] = OMAS_axisTools.updateQuaternion(obj.VIRTUAL.quaternion,globalRates,dt);
            R_new = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = R_new*localNEDState(velocityIndices);
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*globalVelocity_k_plus; % ASSUME INSTANTANEOUS
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.priorState = obj.localState;
            obj.localState = localNEDState;                        % Reassign the obj.localstate
        end
        % UPDATE THE GLOBAL PROPERTIES ( LOCAL AXIS REMAINS STATIC )
        function [obj] = updateGlobalProperties_fixedFrame(obj,dt,localENUState)
           % This function computes the new global parameters for the given
           % agent based on its new state.
            
            % CONSTANT PARAMETERS
            velocityIndices = 7:9;
           
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                obj.VIRTUAL.idleStatus = 1;
                localENUState(7:12) = zeros(6,1); % Freeze the agent
            end   
            
            % /////////////// DEFINE PREVIOUS PARAMETERS //////////////////
            quaternion_k     = obj.VIRTUAL.quaternion;
            velocity_ENU_k_plus = localENUState(velocityIndices,1); 
                        
            % ///////////////// [WORKING: NON-ROTATED] ////////////////////
            quaternion_k_plus = quaternion_k;
            
            %NEW ROTATION MATRIX
            [R_update] = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            localVelocityUpdate = velocity_ENU_k_plus;
            
            % /////////////////// UPDATE GLOBAL POSE //////////////////////
            globalVelocity_k_plus = R_update*localVelocityUpdate;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*globalVelocity_k_plus;
            
            % ///////////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.rotationMatrix = R_update;
            obj.localState = localENUState;  
        end
    end
end

