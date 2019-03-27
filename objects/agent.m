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
        MEMORY;                         % Record of relative scene
        maxSamples = 1;                 % The maximum number of states retained
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
            obj = obj.SetDetectionRadius(inf);
            obj = obj.SetType(OMAS_objectType.agent);
            obj = obj.SetHitBoxType(OMAS_hitBoxType.spherical);
            obj = obj.SetSymbol('diamond');
            obj = obj.SetIdleStatus(logical(false));
            obj = obj.SetRadius(obj.radius);
            
            % Define default sensor model
            [obj.SENSORS] = obj.GetCustomSensorParameters();
            % Initialise memory structure (with 3D varient)
            [obj] = obj.SetBufferSize(obj.maxSamples);
            
            % CHECK FOR USER OVERRIDES
            obj.VIRTUAL = obj.configurationParser(obj.VIRTUAL,varargin); 
            obj = obj.configurationParser(obj,varargin);
        end
        % MAIN 
        function obj = main(obj,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(dt,varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            [headingVector] = obj.GetTargetHeading();
            desiredVelocity = headingVector*obj.nominalSpeed;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
%             [figureHandle] = obj.GetObjectScene(gcf);
            
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [obj] = obj.controller(dt,desiredVelocity);
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    % SET METHODS
    methods 
        % Set the maximum sample number
        function [obj] = SetBufferSize(obj,horizon)
            assert(numel(horizon) == 1 && isnumeric(horizon),'Horizon must be a scalar number of steps.');
            % define the size of the buffer
            obj.maxSamples = horizon;
            % Initialise memory structure (with 3D varient)
            [obj.MEMORY]  = obj.GetMemoryStructure(obj.maxSamples);
        end
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
        % Set the representative radius
        function [obj] = SetRadius(obj,radius)
            % Input sanity check
            assert(isnumeric(radius) && numel(radius) == 1,'Representative radius must be a scalar value.');
            % Set the representative radius
            obj.radius = radius;
            obj.VIRTUAL.radius = radius;        
        end
        % Check dimensionality of object
        function [flag] = is3D(obj)
           flag = obj.VIRTUAL.is3D;
        end
    end
    % ////////////////////// BASIC UPDATE FUNCTIONS ///////////////////////
    methods
        %%%% AGENT UPDATE (IDEAL) %%%%
        % Update from perfect environmental knowledge
        function [obj,obstacleSet,agentSet,waypointSet] = GetAgentUpdate(obj,dt,observedObjects)
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
            
            assert(isnumeric(dt),'First parameter must be a numeric timestep.');
            assert(isstruct(observedObjects),'Second parameter is a vector of observation structures.');
            % Check #1 - Nothing is observed..
            if isempty(observedObjects)
                return
            end
            
            % Update agent memory structure
            for entry = 1:numel(observedObjects)
                % Apply sensor model if there is one
                sensedObject = obj.SensorModel(dt,observedObjects(entry));
                % Update memory structure from measurements
                obj = obj.UpdateMemoryFromObject(sensedObject);
            end
            
            % SORT OBJECT SET BY PRIORITY
            [obj] = obj.SortMemoryByField('priority');
            
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
        % Default sensor model
        function [observedObject] = SensorModel(obj,dt,observedObject)
            
            % Get measurements from a camera
            [d_i,psi_i,theta_i,alpha_i] = obj.GetCameraMeasurements(observedObject);
            % Override ideal parameters with camera model
            observedObject.range = d_i;
            observedObject.heading = psi_i;
            observedObject.elevation = theta_i;
            observedObject.width = alpha_i;
        end
    end
    % ///////////////////////// SENSING FUNCTIONS /////////////////////////
    methods
        % SENSOR MODEL - CAMERA & RANGE FINDER
        function [d_j,psi_j,theta_j,alpha_j] = GetCameraMeasurements(obj,observedObject) 
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
            d_j     = observedObject.range + obj.SENSORS.sigma_rangeFinder*randn(1);     
            psi_j   = observedObject.heading + obj.SENSORS.sigma_camera*randn(1);
            theta_j = observedObject.elevation + obj.SENSORS.sigma_camera*randn(1);
            % OBSERVED RADIUS
            alpha_j = observedObject.width + obj.SENSORS.sigma_camera*randn(1);      % Uncertainty in the angular measurement
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
        function [p] = GetCartesianFromSpherical(r,phi,theta)
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
            r = norm(p);                                        % The range
            phi = atan2(p(2),p(1));                             % The angle made in the azimuth (bearing)
            theta = atan2(p(3),sqrt(p(1).^2 + p(2).^2));        % The elevation (vertical bearing)
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
    %% ////////////////////// SIMULATION INTERFACES ///////////////////////
    methods
        % UPDATE OBJECTIVE
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
        % GET TARGET HEADING VECTOR OF VISIBLE OBJECT
        function [headingVector] = GetTargetHeading(obj,targetObject)
            % This function calculates the heading vector to the current
            % obj.targetWaypoint, or to the provided object with a position
            % field.
            
            % If a target is provided
            if nargin > 1
                targetPosition = obj.GetLastMeasurement(targetObject.objectID,'position');
            elseif ~isempty(obj.targetWaypoint)
                targetPosition = obj.targetWaypoint.position(:,obj.targetWaypoint.sampleNum);
            else
                % Default heading
                if obj.VIRTUAL.is3D
                    targetPosition = [1;0;0];
                else
                    targetPosition = [1;0];
                end
            end
            % Get the vector heading of the target
            headingVector = targetPosition/norm(targetPosition);
        end
        % GET TARGET-GET CONDITION
        function [targetLogical] = GetTargetCondition(obj)
            % Get the current measurements
            currentPosition = obj.targetWaypoint.position(:,obj.targetWaypoint.sampleNum);
            currentRadius   = obj.targetWaypoint.radius(obj.targetWaypoint.sampleNum);
            % Check the current achieve-tolerances on a given target
            targetLogical = 0 > norm(currentPosition) - (currentRadius + obj.VIRTUAL.radius);
        end
    end
    % MEMORY OPERATIONS
    methods
        % //////////////////// MEMORY MANIPULATION ////////////////////////
        % SORT MEMORY   - DEFINED FIELD
        function [obj] = SortMemoryByField(obj,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            assert(isfield(obj.MEMORY,field),'Field must belong to the memory structure.');
            
            % Reorder the memory structure based on fieldname
            [~,ind] = sort([obj.MEMORY.(field)],2,'descend');         % Ordered indices of the object IDs
            % Sort the memory structure 
            obj.MEMORY = obj.MEMORY(ind);
        end
        
        % ///////////////////////// MEMORY I/O ////////////////////////////
        % UPDATE MEMORY - FROM OBSERVATIONS
        function [obj] = UpdateMemoryFromObject(obj,observedObject)
            
            % This program is designed to parse the sensor data recieved into memory
            % for access in main agent code.
            
            % Input sanity check #1
            if numel(observedObject) == 0
                return
            end
                       
            % [TO-DO] REMOVE MEMORY ENTRIES NOT VISIBLE ///////////////////
            % Get the indices of objects still visible
            % visibleIDs = ismember([obj.MEMORY.objectID],[observedObject.objectID]);
            % Remove entries no longer visible
            % obj.MEMORY = obj.MEMORY(visibleIDs);                              % Removes all entries not visible 
            
            % ////////// UPDATE MEMORY WITH NEW ENTRIES AND EXISTING ENTRIES //////////
            % Get the logical indices of ID occurances
            logicalIDIndex = [obj.MEMORY.objectID] == observedObject.objectID;  % Appearance of object ID in memory
            IDoccurances   = sum(logicalIDIndex);
            % Determine memory behaviour
            switch IDoccurances
                case 0 % ///////////////// IS NOT IN MEMORY YET ///////////////////
                    % Override the template if its the first reading
                    if numel(obj.MEMORY) == 1 && obj.MEMORY(1).objectID == 0
                        logicalIDIndex = 1;                                     % The memory address to be overwritten
                    else
                        % Get new memory structure with associated fields
                        obj.MEMORY = vertcat(obj.MEMORY,obj.GetMemoryStructure(obj.maxSamples));        % Append the new structure
                        logicalIDIndex = numel(obj.MEMORY);                     % Address is the end address
                    end
                    % UPDATE MEMORY FIELDS FROM OBSERVATIONS
                    memFields = fieldnames(obj.MEMORY(logicalIDIndex));         % Update fields by dynamic association
                    for entry = 1:numel(memFields)
                        if ~isfield(observedObject,memFields{entry})            % Check the memory field is observable
                            continue                                            % Field cannot be updated from observations
                        end
                        % Data is available
                        if isa(obj.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            % Update circular-buffers
                            data = obj.MEMORY(logicalIDIndex).(memFields{entry});
                            data(:,obj.MEMORY(logicalIDIndex).sampleNum) = observedObject.(memFields{entry});
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                        else
                            % Direct value override
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = observedObject.(memFields{entry});
                        end
                    end
                case 1 % ////////// IS IN MEMORY ALREADY AND SINGUAR //////////////
                    % Indicate new sample
                    obj.MEMORY(logicalIDIndex).sampleNum = obj.MEMORY(logicalIDIndex).sampleNum + 1;
                    % Get the new sample number
                    t = obj.MEMORY(logicalIDIndex).sampleNum;
                    % Update dynamic fields by association
                    memFields = fieldnames(obj.MEMORY(logicalIDIndex));
                    for entry = 1:numel(memFields)
                        % Isolate circular-buffers
                        if isa(obj.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            data = obj.MEMORY(logicalIDIndex).(memFields{entry});
                            data(:,t) = observedObject.(memFields{entry});
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                        end
                    end
                otherwise
                    error('[ERROR] Agent memory structure is distorted.');
            end
            % /////////// DEDUCTIONS FROM THE MEASUREMENTS ////////////
            % [TO-DO] non-proximity based priority for now
            obj.MEMORY(logicalIDIndex).priority = 1/norm(obj.MEMORY(logicalIDIndex).position(:,obj.MEMORY(logicalIDIndex).sampleNum));
        end
        % FROM MEMORY   - WHOLE STRUCTURE BY ID
        function [memStruct] = GetObjectMemoryStruct(obj,objectID)
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            % The memory structure
            memStruct = obj.MEMORY(memoryLogicals);
        end
        % FROM MEMORY   - TRAJECTORY BY FIELD
        function [trajectoryData] = GetTrajectoryByObjectID(obj,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            assert(ichar(field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            trajectoryData = obj.MEMORY(memoryLogicals).(field);
            % Reorder the fields in the data
            trajectoryData = NaN(size(trajectoryData));
            for step = 0:(obj.maxSamples - 1)
                bufferIndex = obj.MEMORY(memoryLogicals).sampleNum - step; % Index in the circular buffer
                trajectoryData(:,ind) = trajectoryData(:,bufferIndex);                   % Write the trajectory in (t,t-1,t-2) order
            end
        end
        % FROM MEMORY   - LAST MEASUREMENT BY FIELD
        function [latestData] = GetLastMeasurementByObjectID(obj,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            assert(ischar(field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            latestData = obj.MEMORY(memoryLogicals).(field);
            % Get the latest samples data from buffer  
            if isa(latestData,'circularBuffer')
                latestData = latestData(:,obj.MEMORY(memoryLogicals).sampleNum);
            else
                latestData = latestData(:,1);
            end
        end
    end
    methods (Static)
        % //////////////////// MEMORY INITIALISATION //////////////////////
        % FROM MEMORY - LAST MEASUREMENT FROM STRUCTURE
        function [latestData] = GetLastMeasurementFromStruct(memStruct,field)
            % Input sanity check
            assert(ischar(field),'Field must be a defined as a string label.');
            
            % Get the associated buffer data
            latestData = memStruct.(field);
            % Get the latest samples data from buffer  
            if isa(latestData,'circularBuffer')
                latestData = latestData(:,memStruct.sampleNum);
            else
                latestData = latestData(:,1);
            end
        end
        % GET EMPTY MEMORY STRUCTURE (3D trajectories)
        function [memStruct]  = GetMemoryStructure(horizonSteps)
            % This function contains a basic agent-memory structure. This
            % is used to retain information on observed objects and maintain 
            % a regular structure.
            
            % Input sanity check
            if nargin < 1
                horizonSteps = 10; % Duration retained in memory
            end
            
            % The fields of memory structure define the fields of the
            % simulation 'observation' structure that are retained.
            
            % Create empty memory structure
            memStruct = struct(...
                'name','',...
                'objectID',uint8(0),...
                'type',OMAS_objectType.misc,...
                'sampleNum',uint8(1),...
                't_lastSeen',[],...
                'position',circularBuffer(NaN(3,horizonSteps)),...
                'velocity',circularBuffer(NaN(3,horizonSteps)),...
                'radius',circularBuffer(NaN(1,horizonSteps)),...
                'range',circularBuffer(NaN(1,horizonSteps)),...
                'heading',circularBuffer(NaN(1,horizonSteps)),...
                'elevation',circularBuffer(NaN(1,horizonSteps)),...
                'width',circularBuffer(NaN(1,horizonSteps)),...
                'geometry',struct('vertices',[],'faces',[],'normals',[],'centroid',[]),...
                'TTC',[],...
                'priority',[]);
        end
    end

    methods
        % PLOT ALL THE OBSERVABLE OBSTACLES IN THE LOCAL FRAME
        function [figureHandle] = GetObjectScene(obj,figureHandle)
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
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
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
            % ADD ANNOTATION
            annotationText = sprintf(' \t%s [ID-%d]',obj.name,obj.objectID);
            text(0,0,0,char(annotationText));       
            
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
            for item = 1:numel(obj.MEMORY)
                % Get the second objects data
                [p_j] = obj.GetLastMeasurementByObjectID(obj.MEMORY(item).objectID,'position');
                [v_j] = obj.GetLastMeasurementByObjectID(obj.MEMORY(item).objectID,'velocity');
                [r_j] = obj.GetLastMeasurementByObjectID(obj.MEMORY(item).objectID,'radius');
                
                if numel(obj.MEMORY(item).geometry.vertices) < 1
                    % DEFINE GEOMETRY AS A SPHERE
                    if obj.VIRTUAL.is3D
                        [geometry] = OMAS_graphics.defineSphere(p_j,r_j);
                    else
                        [geometry] = OMAS_graphics.defineSphere([p_j;0],r_j);
                    end
                else
                    % DEFINE FROM OWN GEOMETRY
                    geometry = obj.MEMORY(item).geometry;                 % Should be in the current frame
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
                switch obj.MEMORY(item).type
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
                annotationText = sprintf(' \t%s [ID-%d]',obj.MEMORY(item).name,obj.MEMORY(item).objectID);
                if obj.VIRTUAL.is3D
                    text(p_j(1),p_j(2),p_j(3),char(annotationText));       % PLOT AGENT VELOCITY
                    q = quiver3(ax,p_j(1),p_j(2),p_j(3),...
                                   v_j(1),v_j(2),v_j(3),'k');              % The velocity vector
                else
                    text(p_j(1),p_j(2),0,char(annotationText));
                    q = quiver3(ax,p_j(1),p_j(2),0,...
                                   v_j(1),v_j(2),0,'k');                   % The velocity vector
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
end

