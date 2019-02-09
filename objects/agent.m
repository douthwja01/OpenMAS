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
        % DEFAULT BEHAVIOUR
        nominalSpeed = 2;               % Default nominal speed (m/s)
        maxSpeed = 4;                   % Defalt maximumal speed (m/s)
        % WAYPOINTS
        targetWaypoint;                 % The current waypoint target
        achievedWaypoints;              % The agents list of locally achieved waypoints
        % DYNAMIC PARAMETERS
        DYNAMICS;
        % SENSOR MEASUREMENT PARAMETERS
        SENSORS = struct('range',inf);  % Only generic parameter is a visual horizon
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
                        
            % OVERRIDE VIRTUAL DEFINITION
            obj.VIRTUAL.type       = OMAS_objectType.agent;
            obj.VIRTUAL.hitBoxType = OMAS_hitBoxType.spherical;
            obj.VIRTUAL.symbol = 'diamond';
            obj.VIRTUAL.detectionRadius = obj.SENSORS.range;  
            obj.VIRTUAL.idleStatus = logical(false);
            
            % CHECK FOR USER OVERRIDES
            obj.VIRTUAL = obj.configurationParser(obj.VIRTUAL,varargin); 
            obj = obj.configurationParser(obj,varargin);
        end
        % ///////////////////// AGENT MAIN CYCLE //////////////////////////
        function obj = main(obj,TIME,varargin)
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
           
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [obj] = obj.controller(dt,desiredVelocity);
        end
    end
    
    % ////////////////////// (3D) SENSORY FUNCTIONS ///////////////////////
    methods
        % WAYPOINT TARGET VECTOR 
        function [headingVector] = getWaypointHeading(obj)
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                headingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
            else 
                %headingVector = [1;0;0];
                if obj.VIRTUAL.is3D
                    headingVector = [1;0;0];
                else
                    headingVector = [1;0];
                end
            end
        end
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
            
            % INFORMATION DIMENSIONALITY
            if obj.VIRTUAL.is3D
                variableIndices = 1:3; 
            else
                variableIndices = 1:2;
            end
            
            % UPDATE THE AGENT'S UNDERSTANDING OF THE ENIVIRONMENT
            for item = 1:length(observedObjects)                
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(observedObjects(item).priority)
                    objectPriority = observedObjects(item).priority;           % Assign priority to memory
                else
                    objectPriority = 1/observedObjects(item).timeToCollision; 
                end
                
                % MEMORY ITEM CONSTRUCTIONS (TRUE TO CURRENT STEP)
                memoryItem = struct('name',observedObjects(item).name,...
                                'objectID',observedObjects(item).objectID,...
                                    'type',observedObjects(item).type,...
                                  'radius',observedObjects(item).radius,...
                                'geometry',observedObjects(item).geometry,... 
                                'position',observedObjects(item).position(variableIndices,1),...
                                'velocity',observedObjects(item).velocity(variableIndices,1),...
                                     'TTC',observedObjects(item).timeToCollision,...
                                'priority',objectPriority,...
                          'globalPosition',observedObjects(item).globalPosition(variableIndices,1),...
                          'globalVelocity',observedObjects(item).globalVelocity(variableIndices,1)); 
                
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
            
            % INFORMATION DIMENSIONALITY
            if obj.VIRTUAL.is3D
                variableIndices = 1:3; 
            else
                variableIndices = 1:2;
            end
            
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
                measuredPosition = measuredPosition(variableIndices,1);
                
                % CALCULATE THE RADIUS FROM THE ANGULAR WIDTH INTERVAL
                measuredRadius = (sin(measuredAlpha/2)/(1-sin(measuredAlpha/2)))*measuredRange;                       % Calculate the apparent radius
                
                % STATE ESTIMATION FUNCTION 
                [priorPosition,priorVelocity] = obj.getLastStateFromMemory(observedObjects(item).objectID);       % Get prior obstacle knowledge
                
                if isempty(priorPosition)
                    % FIRST SIGHT OF OBSTACLE
                    measuredVelocity = NaN(variableIndices,1);
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
        function [range,azimuth,elevation,angularWidth] = getSensorMeasurements(obj,obstacleData) 
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
            range     = obstacleData.range + obj.SENSORS.sigma_rangeFinder*randn(1);     
            azimuth   = obstacleData.azimuthAngle + obj.SENSORS.sigma_camera*randn(1);
            elevation = obstacleData.inclinationAngle + obj.SENSORS.sigma_camera*randn(1);
            % OBSERVED RADIUS
            angularWidth = obstacleData.angularWidth + obj.SENSORS.sigma_camera*randn(1);      % Uncertainty in the angular measurement
        end
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [position,velocity,radius] = getAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            % GET THE LOCAL ABSOLUTE VARIABLES
            position =  obj.localState(1:3,1) + obj.SENSORS.sigma_position*rand(3,1); % Get the absolute position measurement
            velocity =  obj.localState(7:9,1) + obj.SENSORS.sigma_velocity*rand(3,1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            radius = obj.VIRTUAL.radius;
        end
        % SENSOR MODEL - PERFECT SENSING
        function [obj] = getDefaultSensorParameters(obj)
            % This function is designed to populate the SENSOR field with
            % representative sensor uncertainty.
            % BUILD THE SENSOR FACTOR
            obj.SENSORS = struct(...
                'range',inf,...             % Assume the agent has perfect environmental knowledge (m)
                'sigma_position',0.0,...    % Accurate to within 0.5m
                'sigma_velocity',0.0,...    % Accurate to within 0.1m/s
                'sigma_rangeFinder',0.0,... % Accurate to within 0.1m
                'sigma_camera',0.0,...      % One pixel in a 1080p image
                'sampleFrequency',inf);     % Object has perfect precision
        end
    end
        
    % //////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods    
        % CONTROLLER
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % NUMERICAL SANITY CHECK ONE
            assert(~any(isnan(desiredVelocity)) && isnumeric(desiredVelocity),'Requested velocity vector must be a 3D local vector');
            % NUMERICAL SANITY CHECK TWO
            [unitDirection,desiredSpeed] = obj.nullVelocityCheck(desiredVelocity);   
            
            % APPLY SPEED CONSTRAINT
            if abs(desiredSpeed) > obj.maxSpeed
                desiredSpeed = sign(desiredSpeed)*obj.maxSpeed; 
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dHeading] = obj.getVectorHeadingAngles([1;0;0],unitDirection); % Relative heading angles   
            omega = -dHeading/dt;
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[desiredSpeed;0;0],[0;0;omega]);
            obj.localState(1:6)  = obj.localState(1:6) + dt*dX;
            obj.localState(7:12) = dX;
            
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES /////////////// 
            if obj.VIRTUAL.idleStatus
                obj.localState(7:12) = zeros(6,1);
            end
            
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,obj.localState);
        end
        % SPEED AND HEADING (4D) PID CONTROLLER
        function [d_speed,d_heading,obj] = PIDController(obj,targetVelocity)
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
        % BASIC VELOCITY VECTOR CHECK (2D & 3D)
        function [v_unit,v_mag] = nullVelocityCheck(v)
            v_mag = norm(v);
            if v_mag == 0 
                v_unit = zeros(numel(v),1);
                v_unit(1) = 1;      % [1 0 0] % Default to current local forward direction
            else
                v_unit = v/v_mag;
            end
        end
        % DEFINE THE PITCH AND YAW TO MEET TARGET HEADING (2D & 3D)
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
                lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));  % Get the angle, signed by the direction of its cross product
                % GET THE ELEVATION ANGLE
                theta = atan2(U(3),norm(Uh));
            else
            % HANDLE 2D CASE (LOS)
                % GET THE LINE OF SIGHT ANGLE
                lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));
                % GET THE SUDO ELEVATION ANGLE
                theta = 0;
            end
        end
        % CALCULATE THE NEW STATE ESTIMATE
        function [p,v1] = linearStateEstimation(dt,p0,v0,p)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
                      
           % GET THE POSITION
           dX = (p - p0);
           v1 = dX/dt;  % Defines the average velocity
           
           % NO PREVIOUS VELOCITY RECORDED
           if any(isnan(v0))
              p = p;
              v1 = v1;
              return
           end          
        end
        % GET ESTIMATE OF RELATIVE POSITION & VELOCITY
        function [p]    = cartesianFromSpherical(r,phi,theta)
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
           p = [cos(phi)*cos(theta);...
                                sin(phi)*cos(theta);...
                                             sin(theta)]*r;
        end
        % CONVERT CARTESIAN TO SPHERICAL
        function [r,phi,theta]  = sphericalFromCartesian(p)
            r = norm(p);                                            % The range
            phi = atan2(p(2),p(1));                    % The angle made in the azimuth (bearing)
            theta = atan2(p(3),sqrt(p(1).^2 + p(2).^2));  % The elevation (vertical bearing)
        end
        % VALIDATE THE OBSTACLE
        function [tau]  = validateCollision(p,v)
            % CONFIRM WE KNOW THE HEADING OF THE OBSTACLE
            if any(isnan(v))
                tau = -inf;
                return
            end                 
            % Define the time to closest approach (+ve converging)
            tau = -(dot(p,v)/dot(v,v));
        end
    end  
    
    % ////////// AGENT/DERIVATIVES GLOBAL UPDATE/STATE FUNCTIONS //////////
    methods
        % GLOBAL UPDATE - EULER 6DOF(3DOF) TO NED STATE VECTOR
        function [obj] = updateGlobalProperties_NED(obj,dt,eulerState)
            
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
                eulerIndices = 4:6;
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
                eulerIndices = 3;
            else
                error('State notation not recognised');
            end
            
            % EQUIVALENT RATES
            velocity_k_plus   = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt;
            eulerRates_k_plus = (eulerState(eulerIndices) - obj.VIRTUAL.priorState(eulerIndices))/dt;  
                                         
            % NED ROTATION CONVENTION
            % 3D - Rotation order (Z Y X) 
            % 2D - Rotation order (Z)
                            
            % Formulate the conversion
            globalAxisRates = eulerRates_k_plus;
            if numel(globalAxisRates) ~= 3
               globalAxisRates = zeros(3,1) + [-1;0;0]*eulerRates_k_plus;
               velocity_k_plus = [velocity_k_plus;0];
            else
                % GET CONSTANT CONVERSION FROM GLOBAL NED TO GLOBAL ENU
                R_NED_ENU = [1 0 0 ; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
               % Formulate the conversion
               globalAxisRates = W*eulerRates_k_plus;
            end
            
            % INTEGRATE THE GLOBAL QUATERNION 
            [quaternion_k_plus] = OMAS_geometry.integrateQuaternion(obj.VIRTUAL.quaternion,globalAxisRates,dt);
            % NEW ROTATION MATRIX (G>B)            
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(quaternion_k_plus);

            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end    
        % GLOBAL UPDATE - EULER 6DOF(3DOF) TO NED STATE VECTOR ( LOCAL AXIS REMAINS STATIC )
        function [obj] = updateGlobalProperties_NED_fixed(obj,dt,eulerState)
            
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
            else
                error('State notation not recognised');
            end
            
            % EQUIVALENT RATES
            velocity_k_plus   = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt; 
                        
            % NED ROTATION CONVENTION
            % 3D - Rotation order (Z Y X) 
            % 2D - Rotation order (Z)
            % ROTATION RATES ABOUT THE GLOBAL AXES
            
            % /////////////// DEFINE PREVIOUS PARAMETERS //////////////////
            velocity_k_plus = eulerState(velocityIndices,1);    % The velocity in the local NED frame
            
            % Formulate the conversion
            if numel(velocity_k_plus) ~= 3
               velocity_k_plus = [velocity_k_plus;0];
            end
            
            % GET LOCAL NED TO GLOBAL NED ROTATION MATRIX
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(obj.VIRTUAL.quaternion);
            % GET CONSTANT CONVERSION FROM GLOBAL NED TO GLOBAL ENU
            R_NED_ENU = [1 0 0 ; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];

            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_NED_ENU*(R_k_plus*velocity_k_plus);
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end    
    end
    
    % //////////////////// SIMULATION & CORE INTERFACES ///////////////////
    methods
        % PLOT ALL THE OBSERVABLE OBSTACLES IN THE LOCAL FRAME
        function [figureHandle] = getObjectScene(obj,objectSet,figureHandle)
            % FIGURE PREPERATION
            if ~exist('figureHandle','var')
               figureHandle = figure(); 
            else
               cla reset
            end
                        
            ax = get(figureHandle,'CurrentAxes');
            set(ax,'NextPlot','replacechildren')
            hold on; grid on;
            axis equal;
            xlabel('X (m)');
            ylabel('Y (m)');
            zlabel('Z (m)');
            view([-120 35]);
            title(ax,'The agents perspective of its neighbourhood');
            
            % PLOT THE AGENT FOR PERSPECTIVE
            if size(obj.GEOMETRY.vertices,1) < 1
                geometry = OMAS_geometry.defineSphere(zeros(3,1),obj.VIRTUAL.radius);
                % REPRESENT GEOMETRY AS A PATCH
                entityHandle = patch(ax,...
                        'Vertices',geometry.vertices,...
                        'Faces',geometry.faces,...
                        'FaceColor',[0.4940, 0.1840, 0.5560]);
            else
                % PLOT THE AGENTS VERTEX DATA
                entityHandle = patch(ax,...
                        'Vertices',obj.GEOMETRY.vertices,...
                        'Faces',obj.GEOMETRY.faces,...
                        'FaceColor',obj.VIRTUAL.colour);
            end
            % SET THE ENTITY DATA
            set(entityHandle,...
                'EdgeColor','k',...
                'EdgeAlpha',0.3,...
                'FaceLighting','gouraud',...
                'FaceAlpha',0.8,...
                'LineWidth',0.1,...
                'EdgeAlpha',0.7);
            
            % PLOT AGENT VELOCITY
            q = quiver3(0,0,0,1,0,0,'b');
            q.AutoScaleFactor = 1;
            % MOVE THROUGH THE OBSTACLE SET AND PLOT THE OBSTACLE POSITIONS
            for item = 1:length(objectSet)
                entity = objectSet(item);
                if numel(entity.geometry.vertices) < 1
                    % DEFINE GEOMETRY AS A SPHERE
                    if obj.VIRTUAL.is3D
                        [geometry] = OMAS_graphics.defineSphere(entity.position,entity.radius);
                    else
                        [geometry] = OMAS_graphics.defineSphere([entity.position;0],entity.radius);
                    end
                else
                    % DEFINE FROM OWN GEOMETRY
                    geometry = entity.geometry;                 % Should be in the current frame
                end
                % REPRESENT GEOMETRY AS A PATCH
                entityHandle = patch(ax,...
                        'Vertices',geometry.vertices,...
                        'Faces',geometry.faces,...
                        'EdgeColor','k',...
                        'EdgeAlpha',0.2,...
                        'FaceLighting','gouraud',...
                        'FaceAlpha',0.2,...
                        'LineWidth',1);
                % PLOT REPRESENTATION
                switch entity.type
                    case OMAS_objectType.agent
                        set(entityHandle,'FaceColor','b');
                    case OMAS_objectType.obstacle
                        set(entityHandle,'FaceColor','r');
                    case OMAS_objectType.waypoint
                        set(entityHandle,'FaceColor','g');
                    otherwise
                        set(entityHandle,'FaceColor','m');
                end
                % ADD ANNOTATION
                annotationText = sprintf(' \t%s [ID-%s]',entity.name,num2str(entity.objectID));
                if obj.VIRTUAL.is3D
                    text(entity.position(1),entity.position(2),entity.position(3),char(annotationText));                % PLOT AGENT VELOCITY
                    q = quiver3(ax,entity.position(1),entity.position(2),entity.position(3),...
                                   entity.velocity(1),entity.velocity(2),entity.velocity(3),'k'); % The velocity vector
                else
                    text(entity.position(1),entity.position(2),0,char(annotationText));
                    q = quiver3(ax,entity.position(1),entity.position(2),0,...
                                   entity.velocity(1),entity.velocity(2),0,'k'); % The velocity vector
                end
                q.AutoScaleFactor = 1;
            end
            drawnow;
            hold off;
        end
        % APPEND AN FIGURE FRAME TO ANIMATION
        function [figureHandle] = getAnimationFrame(obj,ENV,figureHandle,fileName)
            % INPUT HANDLING
            if nargin < 4
                fileName = sprintf('bodyAxes - %s [ID-%.0f].gif',obj.name,obj.objectID);
            end
            
            % FIGURE PREPERATION
            filePath = strcat(ENV.outputPath,fileName);
            drawnow;                                                       % Ensure all items are plotted before taking the image
            % CAPTURE FIGURE AS FRAME
            frame = getframe(figureHandle);         
            im = frame2im(frame);                                          % Capture the figure as an image
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if ENV.currentStep == 1
                imwrite(imind,cm,filePath,'gif',...
                        'Loopcount',inf,...
                        'DelayTime',ENV.dt,...
                        'Location',[0 0],...
                        'Comment',sprintf('%s[ID-%.0f]',obj.name,obj.objectID)); 
            else
                imwrite(imind,cm,filePath,'gif',...
                        'WriteMode','append',...
                        'DelayTime',ENV.dt,...
                        'Location',[0 0]); 
            end 
            pause(0.1);                                                    % Pause for stability
        end
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
            waypointIndex  = linspace(1,length(priorityVector),length(priorityVector)); % Memory of their position in the waypoint set
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
                        obj.VIRTUAL.idleStatus = logical(true);
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
            waypointGetCondition = 0 > norm(obj.targetWaypoint.position) - (obj.targetWaypoint.radius + obj.VIRTUAL.radius);

            % IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
            if waypointGetCondition && ~any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID))
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
        function [lastMemoryItem]    = getLastEntryFromMemory(obj,objectID)
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
        % AGENT DATA STORAGE
        function [obj] = writeAgentData(obj,TIME,loopIndicator,loopDuration)
            % This function writes a regular DATA structure to the agent
            % .DATA structure to allow it to be easily interpretted.
                        
            % THE AGENT COMPUTATION TIMES
            obj.DATA.indicator(TIME.currentStep) = loopIndicator;
            obj.DATA.dt(TIME.currentStep) = loopDuration;
            obj.DATA.steps(TIME.currentStep) = TIME.currentStep;
            obj.DATA.time(TIME.currentStep) = TIME.currentTime;
        end 
    end
end

