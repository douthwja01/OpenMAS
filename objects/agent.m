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
        MEMORY;                         % Record of last known [objectID;positions;velocities]
        % DEFAULT BEHAVIOUR
        nominalSpeed = 2;               % Default nominal speed (m/s)
        maxSpeed = 4;                   % Default maximumal speed (m/s)
        radius = 0.5;                   % Default radius parameter (m)
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
    
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % CONSTRUCTOR
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
                          
            % Simulation setup
            obj = obj.SetDetectionRadius(obj.SENSORS.range);
            obj = obj.SetType(OMAS_objectType.agent);
            obj = obj.SetHitBoxType(OMAS_hitBoxType.spherical);
            obj = obj.SetSymbol('diamond');
            obj = obj.SetIdleStatus(logical(false));
            
            % CHECK FOR USER OVERRIDES
            obj.VIRTUAL = obj.configurationParser(obj.VIRTUAL,varargin); 
            obj = obj.configurationParser(obj,varargin);
        end
        % MAIN 
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
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            [headingVector] = obj.GetTargetHeading();
            desiredVelocity = headingVector*obj.nominalSpeed;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [obj] = obj.controller(dt,desiredVelocity);
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    % SET METHODS
    methods 
        % Set the idle status
        function [obj] = SetIdleStatus(obj,idleStatus)
            assert(islogical(idleStatus),'The logical status must be a logical.');
            obj.VIRTUAL.idleStatus = idleStatus;
        end
        % Set the maximum speed
        function [obj] = SetNominalSpeed(obj,speed)
           obj.nominalSpeed = speed; 
        end
        % Set the detection radius
        function [obj] = SetDetectionRadius(obj,radius)
            % Input sanity check
            assert(isnumeric(radius) && numel(radius) == 1,'Detection radius must be a scalar value.');
            % Set the detection radius
            obj.VIRTUAL.detectionRadius = radius;
            % Override the sensor parameter to reflect the change
            if isstruct(obj.SENSORS)
                obj.SENSORS.range = radius;  % Set the sensor range value
            end
        end
    end
    % ////////////////////// BASIC UPDATE FUNCTIONS ///////////////////////
    methods
        %%%% AGENT UPDATE (IDEAL) %%%%
        % Update from perfect environmental knowledge
        function [obj,obstacleSet,agentSet,waypointSet] = GetAgentUpdate(obj,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = []; obstacleSet = []; waypointSet = [];             % No observed objects 
            
            % INPUT HANDLING
            if isempty(observedObjects)
                return
            end
            
            % UPDATE THE AGENT'S UNDERSTANDING OF THE ENIVIRONMENT
            for item = 1:length(observedObjects)               
                % UPDATE FROM VIRTUAL SENSORS
                [memoryEntry] = obj.GetObjectUpdateNoSensor(observedObjects(item));
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.UpdateMemoryByEntry(memoryEntry);
            end
            
%             obj = obj.UpdateMemoryFromObjects(observedObjects);
            
            
            % SORT OBJECT SET BY PRIORITY
            [obj.MEMORY] = obj.SortMemoryByField('priority');
            
            % DECERN OBJECT TYPES
            % obj.MEMORY now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet    = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.agent);
            waypointSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.waypoint);
            obstacleSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.obstacle); % Differenciate between the different object types
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.UpdateTargetWaypoint(waypointSet);                    % Update the target waypoint and heading vector
        end
        % GET OBJECT DATA STRUCTURE (IDEAL)
        function [memoryEntry] = GetObjectUpdateNoSensor(obj,observedObject)
            
            % INFORMATION DIMENSIONALITY
            if obj.VIRTUAL.is3D
                variableIndices = 1:3;
            else
                variableIndices = 1:2;
            end
            
            % DEFINE ITS APPARENT PRIORITY
            if ~isnan(observedObject.priority)
                objectPriority = observedObject.priority;           % Assign priority to memory
            else
                objectPriority = 1/observedObject.timeToCollision;
            end
                        
            % MEMORY ITEM CONSTRUCTIONS (TRUE TO CURRENT STEP)
            [memoryEntry] = obj.GetMemoryTemplate_ideal();
            % Allocate the object parmaeters
            memoryEntry.name = observedObject.name;
            memoryEntry.objectID = observedObject.objectID;
            memoryEntry.type     = observedObject.type;
            memoryEntry.globalPosition = observedObject.globalPosition(variableIndices,1);
            memoryEntry.globalVelocity = observedObject.globalVelocity(variableIndices,1);
            % Measured parameters
            memoryEntry.geometry = observedObject.geometry;
            memoryEntry.position = observedObject.position(variableIndices,1);
            memoryEntry.velocity = observedObject.velocity(variableIndices,1);
            memoryEntry.radius   = observedObject.radius;
            memoryEntry.TTC      = observedObject.timeToCollision;
            memoryEntry.priority = objectPriority;
        end        
        
        %%%% AGENT UPDATE (FROM SENSORS) %%%%
        % Update from sensors defined in the 'SENSORS' (USING SPHERICAL DATA & SENSOR MODEL)
        function [obj,obstacleSet,agentSet,waypointSet] = GetSensorUpdate(obj,dt,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % dt              - The unit timstep
            % observationSet  - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = []; obstacleSet = []; waypointSet = [];             % No observed objects
            
            % INPUT HANDLING
            if isempty(observedObjects)
                return
            end
            
            % UPDATE THE AGENT'S UNDERSTANDING OF THE ENIVIRONMENT
            for item = 1:length(observedObjects)
                % Get the update
                [memoryEntry] = obj.GetObjectUpdateFromSensors(dt,observedObjects(item));
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.UpdateMemoryByEntry(memoryEntry);
            end
            
            % SORT OBJECT SET BY PRIORITY
            [obj.MEMORY] = obj.SortMemoryByField('priority');
            
            % DECERN OBJECT TYPES
            % obj.MEMORY now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            
            % PULL MEMORY ITEMS FOR LATER MANIPULAT
            agentSet    = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.agent);
            waypointSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.waypoint);
            obstacleSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.obstacle); % Differenciate between the different object types
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES
            [obj,~] = obj.UpdateTargetWaypoint(waypointSet);  % Update the target waypoint and heading vector
        end
        % GET OBJECT DATA STRUCTURE (SENSORS)
        function [memoryEntry] = GetObjectUpdateFromSensors(obj,dt,observedObject)
            
            % INFORMATION DIMENSIONALITY
            if obj.VIRTUAL.is3D
                variableIndices = 1:3;
            else
                variableIndices = 1:2;
            end
            
            % DEFINE THE SENSED VARIABLES
            % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
            [measuredRange,measuredAzimuth,measuredElevation,measuredAlpha] = obj.GetSensorMeasurments(observedObjects(item));
            % Get the cartesian measurements from the spherical
            measuredPosition = obj.GetCartesianFromSpherical(measuredRange,measuredAzimuth,measuredElevation);  % Calculate new relative position
            % Calculate the radius from the angular measurements
            measuredRadius   = obj.GetObjectRadius(measuredRange,measuredAlpha);   % Calculate the apparent radius
            % Map to 2D/3D
            measuredPosition = measuredPosition(variableIndices,1);
            
            % Pull priori measurements for the given object
            [priorPosition,priorVelocity] = obj.GetLastStateFromMemory(observedObjects(item).objectID);         % Get prior obstacle knowledge
            
            if isempty(priorPosition)
                % FIRST SIGHT OF OBSTACLE
                measuredVelocity = NaN(max(variableIndices),1);
            else
                % CALCULATE THE STATE UPDATE
                [measuredPosition,measuredVelocity] = obj.linearStateEstimation(...
                    dt,priorPosition,priorVelocity,measuredPosition);
            end
            
            % DEFINE ITS APPARENT PRIORITY
            if ~isnan(observedObject.priority)
                objectPriority = observedObject.priority;              % Assign priority to memory
            else
                objectPriority = 1/observedObject.timeToCollision;
            end
            
            % MEMORY ITEM CONSTRUCTIONS (TRUE TO CURRENT STEP)
            [memoryEntry] = obj.GetMemoryTemplate_sensor();
            % Allocate the object parmaeters
            memoryEntry.name = observedObject.name;
            memoryEntry.objectID = observedObject.objectID;
            memoryEntry.type     = observedObject.type;
            memoryEntry.globalPosition = observedObject.globalPosition(variableIndices,1);
            memoryEntry.globalVelocity = observedObject.globalVelocity(variableIndices,1);
            % Measured parameters
            memoryEntry.position = measuredPosition;
            memoryEntry.velocity = measuredVelocity;
            memoryEntry.radius   = measuredRadius;
            memoryEntry.range    = measuredRange;
            memoryEntry.azimuth  = measuredAzimuth;
            memoryEntry.elevation = measuredTheta;
            memoryEntry.alpha    = measuredAlpha;
            memoryEntry.TTC      = observedObjects.timeToCollision;
            memoryEntry.priority = objectPriority;
        end
    end
    % ///////////////////////// SENSING FUNCTIONS /////////////////////////
    methods

        % SENSOR MODEL - CAMERA & RANGE FINDER
        function [d_i,psi_i,theta_i,alpha_i] = GetSensorMeasurements(obj,obstacleData) 
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
            d_i     = obstacleData.range + obj.SENSORS.sigma_rangeFinder*randn(1);     
            psi_i   = obstacleData.azimuthAngle + obj.SENSORS.sigma_camera*randn(1);
            theta_i = obstacleData.inclinationAngle + obj.SENSORS.sigma_camera*randn(1);
            % OBSERVED RADIUS
            alpha_i = obstacleData.angularWidth + obj.SENSORS.sigma_camera*randn(1);      % Uncertainty in the angular measurement
        end
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [p_i,v_i,r_i] = GetAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            % GET THE LOCAL ABSOLUTE VARIABLES
            p_i =  obj.localState(1:3,1) + obj.SENSORS.sigma_position*rand(3,1); % Get the absolute position measurement
            v_i =  obj.localState(7:9,1) + obj.SENSORS.sigma_velocity*rand(3,1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            r_i = obj.VIRTUAL.radius;
        end
    end
    methods (Static)
        % SENSOR MODEL - REPRESENTATIVE SENSING
        function [SENSORS] = GetCustomSensorParameters()
            % This function is designed to populate the SENSOR structure
            % with representative sensing parameters for the assembly of
            % the sensing intervals.
            % BUILD THE SENSOR FACTOR
            SENSORS = struct(...
                'range',inf,...             % Assume the agent has perfect environmental knowledge (m)
                'sigma_position',0.5,...    % Accurate to within 0.5m
                'sigma_velocity',0.1,...    % Accurate to within 0.1m/s
                'sigma_rangeFinder',0.1,... % Accurate to within 0.1m
                'sigma_camera',5.208E-5,... % One pixel in a 1080p image
                'sampleFrequency',inf);     % Object has perfect precision
        end
        % SENSOR MODEL - PERFECT SENSING
        function [SENSORS] = GetDefaultSensorParameters()
            % This function is designed to populate the SENSOR field with
            % perfect sensor capabilities.
            % BUILD THE SENSOR FACTOR
            SENSORS = struct(...
                'range',inf,...             % Assume the agent has perfect environmental knowledge (m)
                'sigma_position',0.0,...    % Perfect position measurement
                'sigma_velocity',0.0,...    % Perfect velocity measurement
                'sigma_rangeFinder',0.0,... % Perfect range acquistion
                'sigma_camera',0.0,...      % Infinte resolution
                'sampleFrequency',inf);     % Object has perfect precision
        end 
    end
    % //////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods         
        % CONTROLLER
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % Input sanity check #1 - Is feasible
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector is not a 2D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.'); 
            
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES ///////////////             
            % Input sanity check #2 - Zero vector
            [unitDirection,desiredSpeed] = obj.nullVelocityCheck(desiredVelocity);    
                     
            % APPLY SPEED CONSTRAINT
            if abs(desiredSpeed) > obj.maxSpeed
                desiredSpeed = sign(desiredSpeed)*obj.maxSpeed; 
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dPsi,dTheta] = obj.GetVectorHeadingAngles([1;0;0],unitDirection); % Relative heading angles   
            dHeading = [0;dTheta;-dPsi];
            omega = dHeading/dt;
            
            % OMIT TRAJECTORY CHANGES IF IDLE
            if obj.VIRTUAL.idleStatus
                omega = zeros(3,1);
                desiredSpeed = 0;
                obj.nominalSpeed = 0;
            end
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[desiredSpeed;0;0],omega);
            obj.localState(1:6)  = obj.localState(1:6) + dt*dX;
            obj.localState(7:12) = dX;
                        
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,obj.localState);
        end
        % SPEED AND HEADING (4D) PID CONTROLLER
        function [obj] = controller_PID(obj,dt,desiredVelocity)
            % This function is desiged to generate feedack from a desired
            % local velocity vector. 
            % INPUTS:
            % targetHeading - The unit heading vector
            % targetSpeed   - The target speed in that direction
            % OUTPUTS:
            % heading_fb    - The feedback on the agent heading
            % speed_fb      - The feedback on the agent speed
            % obj           - The error-updated agent object
            
            % Input sanity check
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector must be a 3D local vector.');
            assert(~any(isnan(desiredVelocity)),'The requested velocity contains NaNs.');
            
            % NUMERICAL SANITY CHECK TWO
            [unitDirection,desiredSpeed] = obj.nullVelocityCheck(desiredVelocity);   
            
            % APPLY SPEED CONSTRAINT
            if abs(desiredSpeed) > obj.maxSpeed
                desiredSpeed = sign(desiredSpeed)*obj.maxSpeed; 
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dPsi,dTheta] = obj.GetVectorHeadingAngles([1;0;0],unitDirection); % Relative heading angles   
            dHeading = [0;dTheta;-dPsi];
            
            % RELATIVE SPEED
            e_speed = desiredSpeed - norm(obj.localState(7:9));            % The speed error 
            controlError = [e_speed;dHeading];                             % -ve in the NED control frame of reference
            
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
            speedFeedback = control_fb(1);                                 % Absolute speed input
            headingFeedback = control_fb(2:4);
            omega = headingFeedback/dt;
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[speedFeedback;0;0],omega);
            obj.localState(1:6)  = obj.localState(1:6) + dt*dX;
            obj.localState(7:12) = dX;
            
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES /////////////// 
            if obj.VIRTUAL.idleStatus
                obj.localState(7:12) = zeros(6,1);
            end
            
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,obj.localState);
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
        function [obj] = GetDynamicParameters(obj)
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
        % CALCULATE THE NEW STATE ESTIMATE
        function [position,velocity] = linearStateEstimation(dt,p0,v0,p1)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
                      
           % GET THE POSITION
           dX = (p1 - p0);
           velocity = dX/dt;    % Defines the average velocity
           position = p1;
           % NO PREVIOUS VELOCITY RECORDED
           if any(isnan(v0))
              return
           end          
        end
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
        function [lambda,theta] = GetVectorHeadingAngles(V,U)
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
        % CALCULATE THE RADIUS 
        function [r] = GetObjectRadius(d,alpha)
        	% Calculate the radius of the object
        	r = (sin(alpha/2)/(1-sin(alpha/2)))*d;     
        end
        % CONVERT SPHERICAL TO CARTESIAN
        function [p]    = GetCartesianFromSpherical(r,phi,theta)
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
        function [r,phi,theta]  = GetSphericalFromCartesian(p)
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
        % GLOBAL UPDATE - DIRECT (AGENT OVERRIDE)
        function [obj] = updateGlobalProperties_direct(obj,p,v,q,X)
            % Under this notation, the state vector already contains the
            % global parameters of the object.
            % INPUTS:
            % globalPosition - 3D global cartesian position given in the ENU coordinate frame.
            % globalVelocity - 3D global cartesian velocity given in the ENU coordinate frame. 
            % quaternion     - The new quaternion pose of body in the global frame.
            % R              - The rotation of the body
            % obj.localState - The previous localState (independant of convention)
            % eulerState     - The new state as reported by the agent            
            
            % Input sanity check
            assert(numel(p) == 3 && size(p,2) == 1,'Global position must be a 3D column vector [3x1].');
            assert(numel(v) == 3 && size(v,2) == 1,'Global velocity must be a 3D column vector [3x1].');
            assert(numel(q) == 4 && size(q,2) == 1,'Global pose must be a 4D quaternion vector [4x1].');
            assert(numel(X) == numel(obj.localState) && size(obj.localState,2) == 1,'The length of the objects state update must match the its local state.');
            
            % Check if the object idle condition is made
            if obj.VIRTUAL.type == OMAS_objectType.agent 
                % Evaluate way-point 
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj = obj.SetIdleStatus(true);
                    obj = obj.SetNominalSpeed(0);
                    v = zeros(3,1);                                        % Freeze the agent
                end   
            end           
            
            % ///////////////// REASSIGN K+1 PARAMETERS //////////////////////////
            % Assign the global parameters
            obj.VIRTUAL.globalPosition = p;                                % Reassign the global position
            obj.VIRTUAL.globalVelocity = v;                                % Reassign the global velocity
            obj.VIRTUAL.quaternion = q;                                    % Reassign the quaternion
            obj.VIRTUAL.R = OMAS_geometry.quaternionToRotationMatrix(q);   % K_plus rotation
            obj.VIRTUAL.priorState = obj.localState;                       % Record the previous state
            obj.localState = X;                                            % Update the current state
        end
    end
    % ///////////////////// SIMULATION I/O INTERFACES /////////////////////
    methods
        % PLOT ALL THE OBSERVABLE OBSTACLES IN THE LOCAL FRAME
        function [figureHandle] = GetObjectScene(obj,objectSet,figureHandle)
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
        function [figureHandle] = GetAnimationFrame(obj,ENV,figureHandle,fileName)
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
    
    %% //////////////////// AGENT MEMORY/WAYPOINT LOGIC ////////////////////
    methods
        % SELECT WAYPOINT
        function [obj,waypointVector] = UpdateTargetWaypoint(obj,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if ~exist('waypointSet','var') || isempty(waypointSet)
                obj.targetWaypoint = [];                % No waypoints are available
                waypointVector = [];                    % Default heading
                obj = obj.SetIdleStatus(logical(true)); % Set idle
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
                        obj = obj.SetIdleStatus(logical(true));            % Set the idle status
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
            waypointGetCondition = obj.GetTargetCondition();
            
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
            [waypointVector] = obj.GetTargetHeading();
        end
        % WAY-POINT CONDITION
        function [targetLogical] = GetTargetCondition(obj)
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            targetLogical = 0 > norm(obj.targetWaypoint.position) - (obj.targetWaypoint.radius + obj.VIRTUAL.radius);
        end
        % GET TARGET HEADING VECTOR 
        function [headingVector] = GetTargetHeading(obj,targetObject)
            % This function calculates the heading vector to the current
            % obj.targetWaypoint, or to the provided object with a position
            % field.
            
            % Input check
            if nargin > 1
                headingVector = targetObject.position/norm(targetObject.position);
            elseif ~isempty(obj.targetWaypoint)
                headingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
            else
                % Default heading
                if obj.VIRTUAL.is3D
                    headingVector = [1;0;0];
                else
                    headingVector = [1;0];
                end
            end
        end 
        
%         % APPEND NEW SENSOR READINGS
%         function [obj] = UpdateMemoryFromObjects(obj,observedObjects)
%             % This function is designed to update the memory structure
%             % based on information of all current objects.
%             % TO DO:
%             % - We need to remove entries for entries no longer visible.
%             % - We need to append state knowledge to known targets..
%             
%             % The agents current list of known objects
%             priorIDSet = [obj.MEMORY(:).objectID];
%             visibleIDSet = [observedObjects(:).objectID];
%             % Reduce the agents knowledge to only those that are seen
%             IDlogicals = ismember(priorIDSet,visibleIDSet);                % true where A are in B
%             % Remove now-unknown targets
%             modifiedIDSet     = priorIDSet(IDlogicals);
%             modifiedMemorySet = obj.MEMORY(IDlogicals);                     % Now contains those that existed and now seen
%             
%             % Get the visible targets not in memory 
%             newVisibleLogicals = ~ismember(visibleIDSet,modifiedIDSet);    % false where A are in B
%             newObjects = observedObjects(newVisibleLogicals);
%             existingObjects = observedObjects(~newVisibleLogicals);
% 
%             % Append data to existing memory structures
%             for item = 1:numel(existingObjects)
%                 % Get the memory structure
%                 [memoryEntry] = obj.GetObjectUpdateNoSensor(existingObjects(item));
%                 
%                 
%             end
%             
%             % Create new memory structures and append
%             for item = 1:numel(newObjects)
%                 [memoryEntry] = obj.GetObjectUpdateNoSensor(newObjects(item));
%                 obj.MEMORY = cat(obj.MEMORY,memoryEntry);
%             end
%             
%             
%         end        
        % UPDATE OBSTACLE KNOWLEDGE
        function [obj] = UpdateMemoryByEntry(obj,newEntry)
            % Updates the agents knowledge based on the new information
            % sensed by the agent. 
            % INPUTS:
            % objectID   - The objects ID
            % objectType - The objects sim type (obstacle or waypoint)
            % radius     - The objects apparent size
            % newState   - The objects estimated state
            % priority   - The objects priority
            
            % Input sanity check
            assert(isstruct(newEntry),'Input must be a memory structure with defined fields.');
            if ~isfield(newEntry,'priority')
                newEntry.priority = 0;
            end

            % CHECK CURRENT KNOWLEDGE STATUS
            if isempty(obj.MEMORY)
                % No knowledge at all
                obj.MEMORY = newEntry; 
                return
            end
            
            % GET THE LOGICAL INDICES OF ID OCCURANCES
            logicalIDIndex = [obj.MEMORY.objectID] == newEntry.objectID;   % Appearance of object ID in memory
            IDoccurances = sum(logicalIDIndex);
            
            % CHECK FOR THE APPEARANCE OF THE OBJECT IN MEMORY
            if IDoccurances > 1
                error('MEMORY DUPLICATE DETECTED.');
            elseif IDoccurances == 0
                % APPEND ENTRY FOR NEW OBJECT
                obj.MEMORY = vertcat(obj.MEMORY,newEntry); % Append to knowledge to other records
            else
                % UPDATE THE MEMORY STRUCTURE WITH THE NEW ENTRY
                obj.MEMORY(logicalIDIndex) = newEntry;                     % Update existing knowledge of  that agent
            end
        end
        % MEMORY SORTER
        function [reorderedMemory]   = SortMemoryByField(obj,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            if ~any(strcmp(fieldnames(obj.MEMORY),field))
                error('Field does not belong to the memory structure.')                
            end
            % Reorder the memory structure based on fieldname
            [~,ind] = sort([obj.MEMORY.(field)],2,'descend');         % Ordered indices of the object IDs
            % Sort the memory structure 
            reorderedMemory = obj.MEMORY(ind);
        end
        % GET OBJECT MEMORY ITEM FOR OBJECT-ID 
        function [lastMemoryItem]    = GetLastEntryFromMemory(obj,objectID)
            % This function is designed to return the complete memory item
            % for a given object ID
            % INPUT:
            % obj.MEMORY - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % lastMemoryItem - The agents last memory record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.MEMORY) || ~any([obj.MEMORY.objectID] == objectID)
                lastMemoryItem = [];
%                 isSuccess = 0;
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.MEMORY.objectID] == objectID;
            lastMemoryItem = obj.MEMORY(IDindex); 
%             isSuccess = 1;
        end
        % GET STATE FOR OBSTACLE ID 
        function [position,velocity] = GetLastStateFromMemory(obj,objectID)
            % This function is designed to return the last known state
            % vector for given obstacle object ID.
            % INPUT:
            % obj.MEMORY - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % lastState - The agents last state record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.MEMORY) || ~any([obj.MEMORY.objectID] == objectID)
                position = [];
                velocity = [];
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.MEMORY.objectID] == objectID;
            position = obj.MEMORY(IDindex).position; 
            velocity = obj.MEMORY(IDindex).velocity; 
        end
    end
    methods (Static)
        % GET EMPTY MEMORY STRUCTURE
        function [entry] = GetMemoryTemplate_ideal()
            % Create emptry memory structure
            entry = struct(...
                'name',[],...
                'objectID',0,...
                'type',[],...
                'position',[],...
                'velocity',[],...
                'radius',[],...
                'range',[],...
                'azimuth',[],...
                'elevation',[],...
                'alpha',[],...
                'TTC',[],...
                'priority',[]);
        end
         % GET EMPTY MEMORY STRUCTURE
        function [entry] = GetMemoryTemplate_sensor()
            % Create emptry memory structure
            entry = struct(...
                'name',[],...
                'objectID',0,...
                'type',[],...
                'position',[],...
                'velocity',[],...
                'radius',[],...
                'range',[],...
                'azimuth',[],...
                'elevation',[],...
                'alpha',[],...
                'TTC',[],...
                'priority',[]);
        end
    end
end

