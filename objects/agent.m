%% THE AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic agent and import this variables 
% into the simulation space for the purpose of multi-vehicle control simulation.
% The agent object is a child of the objectDefintion; the prinicple
% distinctions being:
% sensorRange      - The agent is assumed capable of observing its
%                    surroundings.
% controlFrequency - The frequency at which the control cycle is computed.

% Author: James A. Douthwaite

classdef agent < objectDefinition & agent_tools
%%% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % DEFAULT BEHAVIOUR
        nominalSpeed = 2;               % Default nominal speed (m/s)
        maxSpeed = 4;                   % Default maximumal speed (m/s)
        maxOmega = 0.5;
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
        % Constructor
        function [obj] = agent(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % Call the super class
            obj@objectDefinition(varargin);         

            % //////////////////// SENSOR PARAMETERS //////////////////////
            [obj.SENSORS] = obj.GetDefaultSensorParameters();               % Default sensing
            % /////////////////////////////////////////////////////////////
            
            % Assign defaults
            obj.VIRTUAL  = obj.CreateVIRTUAL();                             % Assign (agent) VIRTUAL structure
            obj.DYNAMICS = obj.CreateDYNAMICS();                            % Assign (agent) VIRTUAL structure
            obj.localState = zeros(12,1);                                   % Assign state
            obj = obj.SetBufferSize(1);                                     % Assign buffersize
            obj = obj.SetDetectionRadius(inf);                              % Assign detection radius        
                        
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            obj = obj.ApplyUserOverrides(varargin);                         % Recursive overrides
            % Re-affirm associated properties
            obj.MEMORY = obj.GetMemoryStructure(obj.maxSamples);            % Initialise memory structure (with 2D varient)             
            obj = obj.SetRadius(obj.radius);                                % Reaffirm radius against .VIRTUAL
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        function [obj] = setup(obj,localXYZvelocity,localXYZrotations)
            
            % Build state vector
            obj.localState = zeros(12,1);
            obj.localState(4:6) = localXYZrotations;
            obj.localState(7:9) = localXYZvelocity;
            
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj = obj.SetVIRTUALparameter('priorState',obj.localState);
        end
        % Main 
        function [obj] = main(obj,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
                       
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                %overHandle = figure('name','testFigure');
                overHandle = gcf;
                ax = gca;
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(ENV,varargin{1});          
            
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
            [obj] = obj.controller(ENV.dt,desiredVelocity);
        end
    end
        
    %% ////////////////////// BASIC UPDATE FUNCTIONS //////////////////////
    methods
        % AGENT UPDATE (IDEAL) - Update from perfect environmental knowledge
        function [obj,obstacleSet,agentSet,waypointSet] = GetAgentUpdate(obj,ENV,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % ENV.dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = []; obstacleSet = []; waypointSet = [];             % No observed objects 
            
            % Check #1 - Nothing is observed..
            if isempty(observedObjects)
                return
            end
                        
            % Input sanity check
            assert(isstruct(ENV),'Expecting environment update structure.');
            assert(isnumeric(ENV.dt),'The time-step must be a numeric timestep.');
            assert(isstruct(observedObjects),'Second parameter is a vector of observation structures.');
            
            % Update agent memory structure
            for entry = 1:numel(observedObjects)
                % Apply sensor model if there is one
                sensedObject = obj.SensorModel(ENV.dt,observedObjects(entry));  
                % Update memory structure from measurements
                obj = obj.UpdateMemoryFromObject(ENV.currentTime,sensedObject);
            end
            
            % SORT OBJECT SET BY PRIORITY
            [obj] = obj.UpdateMemoryOrderByField('priority');
            
            % DECERN OBJECT TYPES
            % obj.MEMORY now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet    = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.agent);
            waypointSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.waypoint);
            obstacleSet = obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.obstacle); % Differenciate between the different object types
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj] = obj.UpdateTargetWaypoint(waypointSet);                    % Update the target waypoint and heading vector
        end
    end
    
    %% //////////////////////// SENSING FUNCTIONS /////////////////////////
    methods
        % SENSOR MODEL - TOP LEVEL (Default)
        function [observedObject] = SensorModel(obj,dt,observedObject)
            % This function provides an overridable method resembling the
            % sensor model through which all objects are processed.
            
            % Input sanity check
            assert(isnumeric(dt) && numel(dt) == 1,'The time step value is invalid.');
            assert(numel(observedObject.position) <= 3,'Expecting a position value to be [2x1].');
            assert(numel(observedObject.velocity) <= 3,'Expecting a velocity value to be [2x1].');
            
            % Get measurements from a camera
            [psi_j,theta_j,alpha_j] = obj.GetCameraMeasurements(observedObject);
            % Override ideal parameters with camera model
            observedObject.heading = psi_j;
            observedObject.elevation = theta_j;
            observedObject.width = alpha_j;
            % Get the range estimate
            [observedObject.range] = obj.GetRangeFinderMeasurements(observedObject);
        end
        % SENSOR MODEL - CAMERA & RANGE FINDER
        function [psi_j,theta_j,alpha_j] = GetCameraMeasurements(obj,observedObject) 
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
            
            % Observed pixel coordinates
            psi_j   = observedObject.heading   + obj.SENSORS.sigma_camera*randn(1);
            theta_j = observedObject.elevation + obj.SENSORS.sigma_camera*randn(1);
            % Observed width
            alpha_j = observedObject.width + obj.SENSORS.sigma_camera*randn(1);      % Uncertainty in the angular measurement
        end
        % SENSOR MODEL - RANGE FINDER
        function [d_j] = GetRangeFinderMeasurements(obj,observedObject)
            % This function takes a simulation data resembling an object
            % and calculates the apparent range to the agent.
            
            % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
            d_j = observedObject.range + obj.SENSORS.sigma_rangeFinder*randn(1);      
        end
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [p_i,v_i,r_i] = GetAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            
            if obj.Is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices  = 1:2;
                velocityIndices = 4:5;
            end
            % GET THE LOCAL ABSOLUTE VARIABLES
            p_i =  obj.localState(positionIndices,1) + obj.SENSORS.sigma_position*rand(numel(positionIndices),1); % Get the absolute position measurement
            v_i =  obj.localState(velocityIndices,1) + obj.SENSORS.sigma_velocity*rand(numel(velocityIndices),1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            r_i = obj.GetRadius();
        end
        % CALCULATE THE NEW STATE ESTIMATE
        function [position,velocity] = linearStateEstimation(obj,dt,p0,v0,p1)
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
    
    %% //////////////////////// DYNAMICS & CONTROL ////////////////////////
    methods         
        % Simple controller
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % Input sanity check #1 - Is feasible
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector is not a 2D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.'); 
            
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES ///////////////             
            % Input sanity check #2 - Zero vector
            [heading,speed] = obj.nullVelocityCheck(desiredVelocity);    
                     
            % Get the relative heading
            [dPsi,dTheta] = obj.GetVectorHeadingAngles([1;0;0],heading); % Relative heading angles
            % The desired heading rate
            omega = ([0;dTheta;-dPsi] - zeros(3,1))/dt;
            
            % Apply kinematic constraints to state changess
            [omega_actual,speed_actual] = obj.ApplyKinematicContraints(dt,omega,speed);
            
            % OMIT TRAJECTORY CHANGES IF IDLE
            if obj.GetIdleStatus()
                omega_actual = zeros(3,1);
                speed_actual = 0;
                obj = obj.SetNominalSpeed(0);
            end
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[speed_actual;0;0],omega_actual);
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
        % Apply the DYNAMICS constraints
        function [omega_actual,speed_actual] = ApplyKinematicContraints(obj,dt,omega,speed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % Determine agent dimension
            if obj.Is3D
                currentVelocity = obj.localState(7:9);
                currentOmega    = obj.localState(10:12);
            else
                currentVelocity = obj.localState(4:5);
                currentOmega    = obj.localState(6); 
            end
            
            % For clarity
            v_lim      = obj.DYNAMICS.maxLinearVelocity;
            dv_lim     = obj.DYNAMICS.maxLinearAcceleration;
            omega_lim  = obj.DYNAMICS.maxAngularVelocity;
            dOmega_lim = obj.DYNAMICS.maxAngularAcceleration;
            
            % /////////////// Constrained heading change //////////////////
            % The implied angular acceleration
            omega_dot = (omega - currentOmega)/dt;
            % Constrain against heading acceleration limits
            bounded_omega_dot = boundValue(omega_dot,-dOmega_lim,dOmega_lim);
            % Calculate actual heading velocity limits
            omega = currentOmega + bounded_omega_dot*dt;
            % Constrain against heading rate limits
            omega_actual = boundValue(omega,-omega_lim,omega_lim);
            
            % ///////////////// Constrained speed change //////////////////
            % The current speed
            currentSpeed = norm(currentVelocity);
            % The proposed acceleration
            s_dot = (speed - currentSpeed)/dt; 
            % Constrain against linear acceleration limits
            bounded_s_dot = boundValue(s_dot,-dv_lim(1),dv_lim(1));
            % Calcullate the actual speed
            speed = currentSpeed + bounded_s_dot*dt;
            % Interesect the velocity with permissible velocities
            speed_actual = boundValue(speed,-v_lim(1),v_lim(1));
        end
        % Create the DYNAMICS structure
        function [DYNAMICS] = CreateDYNAMICS(obj)
            % Define basic DYNAMIC container
            [DYNAMICS] = obj.CreateDYNAMICS_default();
        end
        % Create the unconstrained DYNAMICS
        function [DYNAMICS] = CreateDYNAMICS_default(obj)
          % Define an infinite "throw" range by default
            DYNAMICS = struct();
            if obj.Is3D
                % States defined as [x,y,z,phi,theta,psi,dx,dy,dz,dphi,dtheta,dpsi]^T
                DYNAMICS.maxLinearVelocity      = ones(3,1)*obj.maxSpeed;  % Limits on the agents linear velocity 
                DYNAMICS.maxLinearAcceleration  = inf(3,1);                % Limits on the agents linear acceleration
                DYNAMICS.maxAngularVelocity     = ones(3,1)*obj.maxOmega;  % Limits on the agents angular velocity
                DYNAMICS.maxAngularAcceleration = inf(3,1);                % Limits on the agents angular acceleration
            else
                % States defined as [x,y,psi,dx,dy,dpsi]^T
                DYNAMICS.maxLinearVelocity      = ones(2,1)*obj.maxSpeed;  % Limits on the agents linear velocity 
                DYNAMICS.maxLinearAcceleration  = inf(2,1);                % Limits on the agents linear acceleration
                DYNAMICS.maxAngularVelocity     = obj.maxOmega;            % Limits on the agents angular velocity
                DYNAMICS.maxAngularAcceleration = inf(1,1);                % Limits on the agents angular acceleration
            end
        end
    end
    
    %% ////////////////////// GENERAL STATIC METHODS //////////////////////
    methods (Static)
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
        function [r]    = GetRadiusFromAngularWidth(d,alpha)
        	% Calculate the radius of the object
        	r = (sin(alpha/2)/(1-sin(alpha/2)))*d;     
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
    end 
    
    %% ////////////////// GLOBAL UPDATE/STATE FUNCTIONS ///////////////////
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
            
            % No update required
%             if obj.IsIdle()   % if agents are dynamic, we want to see
%             them settle.
%                 return
%             end
            
            % Input sanity check
            assert(numel(p) == 3 && size(p,2) == 1,'Global position must be a 3D column vector [3x1].');
            assert(numel(v) == 3 && size(v,2) == 1,'Global velocity must be a 3D column vector [3x1].');
            assert(numel(q) == 4 && size(q,2) == 1,'Global pose must be a 4D quaternion vector [4x1].');
            assert(numel(X) == numel(obj.localState) && size(obj.localState,2) == 1,'The length of the objects state update must match the its local state.');
            
            % Check if the object idle condition is made
            if obj.VIRTUAL.type == OMAS_objectType.agent 
                % Evaluate way-point 
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj = obj.SetVIRTUALparameter('idleStatus',true);
%                     v = zeros(3,1);                                        % Freeze the agent
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
    % Only the objectDefinition class has access
    methods (Static, Access = private)  
        % Create agent VIRTUAL structure (Override the objectDefinition)
        function [VIRTUAL]  = CreateVIRTUAL()
            % This is function that assembles the template desciption of an
            % agent type object.
            
            % Define the VIRTUAL structure
            VIRTUAL = struct();
            VIRTUAL.type = OMAS_objectType.agent;
            VIRTUAL.hitBoxType = OMAS_hitBoxType.spherical;
            VIRTUAL.detectionRadius = inf;
            VIRTUAL.radius = 0.5;                      % Diameter of 1m
            VIRTUAL.colour = rand(1,3,'single');       % Virtual colour (for plotting)
            VIRTUAL.symbol = 'diamond';                % Representative symbol
            VIRTUAL.globalPosition = [0;0;0];          % Global Cartesian position
            VIRTUAL.globalVelocity = [0;0;0];          % Global Cartesian velocity
            VIRTUAL.quaternion = [1;0;0;0];            % Global quaternion pose
            VIRTUAL.idleStatus = false;                % Object idle logical
            VIRTUAL.is3D = true;                       % Object operates in 3D logical
            VIRTUAL.priorState = [];                            
        end
    end
    %% ///////////////// VISUALISATION AND GENERAL TOOLS //////////////////
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
                    if obj.Is3D
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
                if obj.Is3D
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
    %% ///////////////////////// GET/SET METHODS //////////////////////////
    methods 
        % Get the idle status
        function [flag] = GetIdleStatus(obj)
            flag = obj.GetVIRTUALparameter('idleStatus');
        end
        % Set the idle status
        function [obj]  = SetIdleStatus(obj,idleStatus)
            % Input sanity check
            assert(islogical(idleStatus),'The logical status must be a logical.');
            obj = obj.SetVIRTUALparameter('idleStatus',idleStatus);
        end
        % Set the maximum speed
        function [obj]  = SetMaxSpeed(obj,v_max)
           assert(numel(v_max) == 1,'Expecting a scalar speed value.'); 
            obj.maxSpeed = v_max;
            obj.DYNAMICS.maxLinearVelocity = v_max; 
        end
        % Set the nominal speed
        function [obj]  = SetNominalSpeed(obj,speed)
           obj.nominalSpeed = speed; 
        end
        % Set the detection radius
        function [obj]  = SetDetectionRadius(obj,radius)
            % Input sanity check
            assert(isnumeric(radius) && numel(radius) == 1,'Detection radius must be a scalar value.');
            % Set the detection radius
            obj = obj.SetVIRTUALparameter('detectionRadius',radius);
            % Override the sensor parameter to reflect the change
            if isstruct(obj.SENSORS)
                % Set the sensor range value
                obj.SENSORS.range = radius;  
            end
        end
        % Get the representative radius
        function [r]    = GetRadius(obj)
            r = obj.GetVIRTUALparameter('radius');
            assert(obj.radius == r,"Radius parameter is different to the VIRTUAL representation.");
        end
        % Set the representative radius
        function [obj]  = SetRadius(obj,radius)
            % Input sanity check
            assert(isnumeric(radius) && numel(radius) == 1,'Representative radius must be a scalar value.');
            % Set the representative radius
            obj.radius = radius;
            obj.SetVIRTUALparameter('radius',radius);
        end
    end
end

