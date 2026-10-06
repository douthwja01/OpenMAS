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
    %% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % DEFAULT BEHAVIOUR
        v_nominal = 2;                  % Default nominal linear speed (m/s)
        v_max  = 4;                     % Default maximal linear speed (m/s)
        w_max  = 0.5;                   % Default maximal angular speed (rad/s)
        detectionRadius = inf;
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
        function [this] = agent(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % this     - The constructed object
            
            % Call the super class
            this@objectDefinition(varargin);
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
            [this.SENSORS] = this.GetDefaultSensorParameters();               % Default sensing
            % /////////////////////////////////////////////////////////////
            this.SetGLOBAL(this.CreateGLOBAL());    % Assign (agent) GLOBAL structure
            this.DYNAMICS = this.CreateDYNAMICS();  % Assign (agent) DYNAMICS structure
            this.SetBufferSize(1);                  % Assign default buffer size 
            
            % Assign defaults
            this.radius = 0.5;                      % Default radius (m)
            this.localState = zeros(12,1);          % Assign state

            % //////////////// Check for user overrides ///////////////////
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            this = this.ApplyUserOverrides(varargin);                         % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup (agent default)
        function [this] = setup(this,localXYZvelocity,localXYZrotations)   % [x y z phi theta psi]
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % INPUTS:
            
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
            
            switch this.Is3D()
                case false 
                    % Initialises the state vector [x y psi ; dx dy dpsi]
                    this.setup_2DVelocities(localXYZvelocity,localXYZrotations);
                case true
                    % Initialises the state vector [x y z phi theta psi ; dx dy dz dphi dtheta dpsi]
                    this.setup_3DVelocities(localXYZvelocity,localXYZrotations)
                otherwise
                    error('Object dimensionality unknown');
            end
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % this      - The updated project
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                %overHandle = figure('name','testFigure');
                overHandle = gcf;
                ax = gca;
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            [headingVector] = this.GetTargetHeading();
            desiredVelocity = headingVector*this.v_nominal;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [this] = this.Controller(ENV.dt,desiredVelocity);
        end
    end
    
    %% ////////////////////// BASIC UPDATE FUNCTIONS //////////////////////
    methods
        % AGENT UPDATE (IDEAL) - Update from perfect environmental knowledge
        function [this,obstacleSet,agentSet,waypointSet] = GetAgentUpdate(this,ENV,observedobjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
            % INPUTS:
            % ENV.dt              - The unit timstep
            % observedobjects - The full observed object structure
            % OUTPUTS:
            % this             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = []; obstacleSet = []; waypointSet = [];             % No observed objects
            
            % Check #1 - Nothing is observed..
            if isempty(observedobjects)
                return
            end
            
            % Input sanity check
            assert(isstruct(ENV),'Expecting environment update structure.');
            assert(isnumeric(ENV.dt),'The time-step must be a numeric timestep.');
            assert(isstruct(observedobjects),'Second parameter is a vector of observation structures.');
            
            % Update agent memory structure
            for entry = 1:numel(observedobjects)
                % Apply sensor model if there is one
                sensedobject = this.SensorModel(ENV.dt,observedobjects(entry));
                % Update memory structure from measurements
                this = this.UpdateMemoryFromObject(ENV.currentTime,sensedobject);
            end
            
            % SORT object SET BY PRIORITY
            [this] = this.UpdateMemoryOrderByField('priority');
            
            % DECERN object TYPES
            % this.MEMORY now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet    = this.MEMORY([this.MEMORY.type] == OMAS_objectType.agent);
            waypointSet = this.MEMORY([this.MEMORY.type] == OMAS_objectType.waypoint);
            obstacleSet = this.MEMORY([this.MEMORY.type] == OMAS_objectType.obstacle); % Differenciate between the different object types
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES
            [this] = this.UpdateTargetWaypoint(waypointSet);                    % Update the target waypoint and heading vector
        end
    end
    
    %% //////////////////////// SENSING FUNCTIONS /////////////////////////
    methods
        % SENSOR MODEL - TOP LEVEL (Default)
        function [observedobject] = SensorModel(this,dt,observedobject)
            % This function provides an overridable method resembling the
            % sensor model through which all objects are processed.
            
            % Input sanity check
            assert(isnumeric(dt) && numel(dt) == 1,'The time step value is invalid.');
            assert(numel(observedobject.position) <= 3,'Expecting a position value to be [2x1].');
            assert(numel(observedobject.velocity) <= 3,'Expecting a velocity value to be [2x1].');
            
            % Get measurements from a camera
            [psi_j,theta_j,alpha_j] = this.GetCameraMeasurements(observedobject);
            % Override ideal parameters with camera model
            observedobject.heading = psi_j;
            observedobject.elevation = theta_j;
            observedobject.width = alpha_j;
            % Get the range estimate
            [observedobject.range] = this.GetRangeFinderMeasurements(observedobject);
        end
        % SENSOR MODEL - CAMERA & RANGE FINDER
        function [psi_j,theta_j,alpha_j] = GetCameraMeasurements(this,observedobject)
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
            psi_j   = observedobject.heading   + this.SENSORS.sigma_camera*randn(1);
            theta_j = observedobject.elevation + this.SENSORS.sigma_camera*randn(1);
            % Observed width
            alpha_j = observedobject.width + this.SENSORS.sigma_camera*randn(1);      % Uncertainty in the angular measurement
        end
        % SENSOR MODEL - RANGE FINDER
        function [d_j] = GetRangeFinderMeasurements(this,observedobject)
            % This function takes a simulation data resembling an object
            % and calculates the apparent range to the agent.
            
            % EMULATE MEASURED POSITIONAL VARIABLES (measured = range(true) + distortion)
            d_j = observedobject.range + this.SENSORS.sigma_rangeFinder*randn(1);
        end
        % SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [p_i,v_i,r_i] = GetAgentMeasurements(this)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in this.localState).
            
            if this.Is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices  = 1:2;
                velocityIndices = 4:5;
            end
            % GET THE LOCAL ABSOLUTE VARIABLES
            p_i =  this.localState(positionIndices,1) + this.SENSORS.sigma_position*rand(numel(positionIndices),1); % Get the absolute position measurement
            v_i =  this.localState(velocityIndices,1) + this.SENSORS.sigma_velocity*rand(numel(velocityIndices),1); % Get the absolute velocity measurement
            % RADIUS IS ASSUMED KNOWN
            r_i = this.radius();
        end
        % CALCULATE THE NEW STATE ESTIMATE
        function [position,velocity] = linearStateEstimation(this,dt,p0,v0,p1)
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
                'sampleFrequency',inf);     % object has perfect precision
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
                'sampleFrequency',inf);     % object has perfect precision
        end
    end
    
    %% //////////////////////// DYNAMICS & CONTROL ////////////////////////
    methods
        % Simple controller
        function [this] = Controller(this,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % Input sanity check #1 - Is feasible
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector is not a 2D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.');
            
            % ///////////// UPDATE object GLOBAL PROPERTIES ///////////////
            % Input sanity check #2 - Zero vector
            [heading,speed] = this.nullVelocityCheck(desiredVelocity);
            
            % Get the relative heading
            [dPsi,dTheta]   = this.GetVectorHeadingAngles([1;0;0],heading); % Relative heading angles
            % The desired heading rate
            omega = ([0;dTheta;-dPsi] - zeros(3,1))/dt;
            
            % Apply kinematic constraints to state changess
            [omega_actual,speed_actual] = this.ApplyKinematicContraints(dt,omega,speed);
            
            % OMIT TRAJECTORY CHANGES IF IDLE
            if this.IsIdle()
                omega_actual = zeros(3,1);
                speed_actual = 0;
                this.v_nominal = 0;
            end
            
            dX = [speed_actual;0;0;omega_actual];
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = this.SimpleDynamics(this.localState(1:6),[speed_actual;0;0],omega_actual);
            this.localState(1:6)  = this.localState(1:6) + dt*dX;
            this.localState(7:12) = dX;
            
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [this] = this.GlobalUpdate_3DVelocities(dt,this.localState);
        end
        % SPEED AND HEADING (4D) PID CONTROLLER
        function [this] = Controller_PID(this,dt,desiredVelocity)
            % This function is desiged to generate feedack from a desired
            % local velocity vector.
            % INPUTS:
            % targetHeading - The unit heading vector
            % targetSpeed   - The target speed in that direction
            % OUTPUTS:
            % heading_fb    - The feedback on the agent heading
            % speed_fb      - The feedback on the agent speed
            % this           - The error-updated agent object
            
            % Input sanity check
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector must be a 3D local vector.');
            assert(~any(isnan(desiredVelocity)),'The requested velocity contains NaNs.');
            
            % NUMERICAL SANITY CHECK TWO
            [unitDirection,desiredSpeed] = this.nullVelocityCheck(desiredVelocity);
            
            % APPLY SPEED CONSTRAINT
            if abs(desiredSpeed) > this.v_max
                desiredSpeed = sign(desiredSpeed)*this.v_max;
            end
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dPsi,dTheta] = this.GetVectorHeadingAngles([1;0;0],unitDirection); % Relative heading angles
            dHeading = [0;dTheta;-dPsi];
            
            % RELATIVE SPEED
            e_speed = desiredSpeed - norm(this.localState(7:9));            % The speed error
            controlError = [e_speed;dHeading];                             % -ve in the NED control frame of reference
            
            % TUNING PARAMETERS
            Kp_linear = 0.8;
            Kd_linear = 0;
            Kp_angular = 1;
            Kd_angular = 0;
            
            % CALCULATE THE CONTROL FEEDBACK
            control_fb = diag([Kp_linear Kp_angular Kp_angular Kp_angular])*controlError + ...
                diag([Kd_linear Kd_angular Kd_angular Kd_angular])*(controlError - this.priorError);
            % REMEMBER PREVIOUS ERROR
            this.priorError = controlError;
            
            % PARSE THE CONTROL FEEDBACK FOR EXTERNAL USE
            speedFeedback = control_fb(1);                                 % Absolute speed input
            headingFeedback = control_fb(2:4);
            omega = headingFeedback/dt;
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = this.SimpleDynamics(this.localState(1:6),[speedFeedback;0;0],omega);
            this.localState(1:6)  = this.localState(1:6) + dt*dX;
            this.localState(7:12) = dX;
            
            % ///////////// UPDATE object GLOBAL PROPERTIES ///////////////
            if this.GLOBAL.idleStatus
                this.localState(7:12) = zeros(6,1);
            end
            
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [this] = this.GlobalUpdate_3DVelocities(dt,this.localState);
        end
        % Apply the DYNAMICS constraints
        function [omega_actual,speed_actual] = ApplyKinematicContraints(this,dt,omega,speed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % Determine agent dimension
            if this.Is3D
                currentVelocity = this.localState(7:9);
                currentOmega    = this.localState(10:12);
            else
                currentVelocity = this.localState(4:5);
                currentOmega    = this.localState(6);
            end
            
            % For clarity
            v_lim      = this.DYNAMICS.maxLinearVelocity;
            dv_lim     = this.DYNAMICS.maxLinearAcceleration;
            omega_lim  = this.DYNAMICS.maxAngularVelocity;
            dOmega_lim = this.DYNAMICS.maxAngularAcceleration;
            
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
        function [DYNAMICS] = CreateDYNAMICS(this)
            % Define basic DYNAMIC container
            [DYNAMICS] = this.CreateDYNAMICS_default();
        end
        % Create the default DYNAMICS constraint structure
        function [DYNAMICS] = CreateDYNAMICS_default(this)
            % Define an infinite "throw" range by default
            DYNAMICS = struct();
            % Define kinematic constraints
            if this.Is3D
                % States defined as [x,y,z,phi,theta,psi,dx,dy,dz,dphi,dtheta,dpsi]^T
                DYNAMICS.maxLinearVelocity      = ones(3,1)*this.v_max;    % Limits on the agents linear velocity
                DYNAMICS.maxLinearAcceleration  = inf(3,1);                % Limits on the agents linear acceleration
                DYNAMICS.maxAngularVelocity     = ones(3,1)*this.w_max;    % Limits on the agents angular velocity
                DYNAMICS.maxAngularAcceleration = inf(3,1);                % Limits on the agents angular acceleration
            else
                % States defined as [x,y,psi,dx,dy,dpsi]^T
                DYNAMICS.maxLinearVelocity      = ones(2,1)*this.v_max;    % Limits on the agents linear velocity
                DYNAMICS.maxLinearAcceleration  = inf(2,1);                % Limits on the agents linear acceleration
                DYNAMICS.maxAngularVelocity     = this.w_max;              % Limits on the agents angular velocity
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
    
    %% ///////////////// VISUALISATION AND GENERAL TOOLS //////////////////
    methods
        % PLOT ALL THE OBSERVABLE OBSTACLES IN THE LOCAL FRAME
        function [figureHandle] = GetobjectScene(this,figureHandle)
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
            if size(this.GEOMETRY.vertices,1) < 1
                geometry = OMAS_geometry.defineSphere(zeros(3,1),this.radius);
                % REPRESENT GEOMETRY AS A PATCH
                entityHandle = patch(ax,...
                    'Vertices',geometry.vertices,...
                    'Faces',geometry.faces,...
                    'FaceColor',[0.4940, 0.1840, 0.5560]);
            else
                % PLOT THE AGENTS VERTEX DATA
                entityHandle = patch(ax,...
                    'Vertices',this.GEOMETRY.vertices,...
                    'Faces',this.GEOMETRY.faces,...
                    'FaceColor',this.GLOBAL.colour);
            end
            % ADD ANNOTATION
            annotationText = sprintf(' \t%s [ID-%d]',this.name,this.objectID);
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
            for item = 1:numel(this.MEMORY)
                % Get the second objects data
                [p_j] = this.GetLastMeasurementByobjectID(this.MEMORY(item).objectID,'position');
                [v_j] = this.GetLastMeasurementByobjectID(this.MEMORY(item).objectID,'velocity');
                [r_j] = this.GetLastMeasurementByobjectID(this.MEMORY(item).objectID,'radius');
                
                if numel(this.MEMORY(item).geometry.vertices) < 1
                    % DEFINE GEOMETRY AS A SPHERE
                    if this.Is3D
                        [geometry] = OMAS_graphics.defineSphere(p_j,r_j);
                    else
                        [geometry] = OMAS_graphics.defineSphere([p_j;0],r_j);
                    end
                else
                    % DEFINE FROM OWN GEOMETRY
                    geometry = this.MEMORY(item).geometry;                 % Should be in the current frame
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
                switch this.MEMORY(item).type
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
                annotationText = sprintf(' \t%s [ID-%d]',this.MEMORY(item).name,this.MEMORY(item).objectID);
                if this.Is3D
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
        function [figureHandle] = GetAnimationFrame(this,ENV,figureHandle,fileName)
            % INPUT HANDLING
            if nargin < 4
                fileName = sprintf('bodyAxes - %s [ID-%.0f].gif',this.name,this.objectID);
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
                    'Comment',sprintf('%s[ID-%.0f]',this.name,this.objectID));
            else
                imwrite(imind,cm,filePath,'gif',...
                    'WriteMode','append',...
                    'DelayTime',ENV.dt,...
                    'Location',[0 0]);
            end
            pause(0.1);                                                    % Pause for stability
        end
        % AGENT DATA STORAGE
        function [this] = writeAgentData(this,TIME,loopIndicator,loopDuration)
            % This function writes a regular DATA structure to the agent
            % .DATA structure to allow it to be easily interpretted.
            
            % THE AGENT COMPUTATION TIMES
            this.DATA.indicator(TIME.currentStep) = loopIndicator;
            this.DATA.dt(TIME.currentStep) = loopDuration;
            this.DATA.steps(TIME.currentStep) = TIME.currentStep;
            this.DATA.time(TIME.currentStep) = TIME.currentTime;
        end
    end
    %% ///////////////////////// GET/SET METHODS //////////////////////////
    methods
        % Set the maximum linear speed
        function set.v_max(this,v)
            assert(numel(v) == 1,'Expecting a scalar linear rate.');
            this.v_max = v;
            this.DYNAMICS.maxLinearVelocity(:) = v;
        end
        % Set the maximum angular speed
        function set.w_max(this,w)
            assert(numel(w) == 1,'Expecting a scalar angular rate.');
            this.w_max = w;
            this.DYNAMICS.maxAngularVelocity(:) = w;
        end
        % Set the detection radius
        function set.detectionRadius(this,r)
            assert(isnumeric(r) && numel(r) == 1,'Detection radius must be a scalar value.');
            % Override the SENSORS parameter to reflect the change
            if isstruct(this.SENSORS) && isfield(this.SENSORS,'range')
                this.SENSORS.range = r;                                    % Set the sensor range value
            end            
            % Assign the value to the GLOBAL structure
            this.SetGLOBAL('detectionRadius',r);                           % Set the detection radius
        end
        % User update (Override)
        function [this] = ApplyUserOverrides(this,pairArray)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called in the class constructor
            
            % Input sanity check #1
            if nargin < 2 || numel(pairArray) == 0
                return
            end
            
            % Call the all parameter overrider
            [this] = GetParameterOverrides_recursive(this,pairArray);
            
            % Preform dependant value checks
            this.MEMORY = this.CreateMEMORY(this.maxSamples);        % Initialise memory structure (with 2D varient)
            
            if this.detectionRadius ~= this.SENSORS.range
                % The detection radius is determined by SENSORS.range
                this.SENSORS.range = this.detectionRadius;
            end
        end
    end
end

