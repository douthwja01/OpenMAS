classdef quadcopter < agent
    properties
        % GLOBAL.position   - The position of the object in the global axes
        % GLOBAL.velocity   - The velocity of the object in the global axes
        % GLOBAL.quaternion - The quaternion representing the earth to rotated body
        
        % DYNAMICS - All the models parameters are held in the DYNAMICS
        %            container field.
    end
    %% ////////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function [this] = quadcopter(varargin)
            % Call the super class
            this@agent(varargin);
            
            % Create the dynamics of a quadcopter
            this.DYNAMICS = this.CreateDYNAMICS();
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a quadrotor
            % object.
            
            % For a state vector defined as x_k:
            % [x,y,z,dx,dy,dz,R11,R12,R13,R21,R22,R23,R31,R32,R33,wx,wy,wz]
            
            % THIS MODEL OPERATES IN THE GLOBAL FRAME
            p0 = this.GetGLOBAL('position');   % True global position
            v0 = this.GetGLOBAL('velocity');   % True global velocity
            % GET THE ROTATION MATRIX fixed local axes -> rotated global pose
            R0 = OMAS_geometry.eulersToRotationMatrix(localXYZrotations);
            % Build an initial (global) state vector
            vecR0   = reshape(R0',9,1);        % Defines local vector to global vector
            omega0  = zeros(3,1);
            x0      = [p0;v0;vecR0;omega0];
            % ASSIGN THE LOCAL FRD STATE
            this.SetGLOBAL('priorState',x0);
            this.localState = x0;
        end
        % Main
        function [this] = main(this,ENV,varargin)
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % Update the agent with the environmental data
            [this,~,~] = this.GetAgentUpdate(ENV,varargin{1});
            % /////////////////////////////////////////////////////////////
            
            % The state reference
            desiredPosition = [0;0;0];
            
            % Call the controller loop
            this = this.Controller(ENV,desiredPosition);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            this.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                '$p$ (rad/s)','$q$ (rad/s)','$r$ (rad/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = [this.localState(4:6);this.localState(16:end)];          % Record the control inputs
        end
    end
    %% //////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % The quadcopter controller
        function [this] = Controller(this,ENV,y_k)
            
            % Input sanity checks
            assert(IsColumn(y_k,3),"Expecting a numeric velocity vector.");
            assert(1 == size(this.localState,2),"The length of the objects state update must match the its local state.");
            
            % ///////////// UPDATE THE LOCAL STATE VECTOR /////////////////
            x_k_plus = this.UpdateLocalState(ENV,this.localState,y_k);
            % Update the global properties
            this = this.GlobalUpdate(x_k_plus);
        end
        % Get the state update (using ode45)
        function [x_k_plus] = UpdateLocalState(this,ENV,x_k,y_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Integrate across the time delta "dt"
            opts = odeset('RelTol',1e-2,'AbsTol',ENV.dt*1E-1);
            
            [~,Xset] = ode45(@(t,X) this.ClosedLoopDynamics_position(X,y_desired),[0 ENV.dt],x_k,opts);
%             [~,Xset] = ode45(@(t,X) this.ClosedLoopDynamics_velocity(X,y_desired),[0 ENV.dt],x_k,opts);
            
            % Assign the integral state
            x_k_plus = Xset(end,:)';
        end
        % Global update from the new state
        function [this] = GlobalUpdate(this,x_k_plus)
            % This function updates the global structure from the new state
            % definition: [p_k;v_k;vecR_k;omega_k]
            
            % Extract the rotation matrix components
            R_k_plus = reshape(x_k_plus(7:15),3,3);
            % Define the new quaternion from R
            q_k_plus = OMAS_geometry.rotationMatrixToQuaternion(R_k_plus');
            
            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            [this] = this.GlobalUpdate_direct(...
                x_k_plus(1:3),...   % The global position is within the state
                x_k_plus(4:6),...   % The global velocity is within the state
                q_k_plus);          % Append the global quaternion
            
            % Ensure the local state is re-assigned
            this.localState = x_k_plus;
        end
    end
    % Dynamic functions
    methods
        % Quadcopter Dynamics (Closed-loop)
        function [dxdt] = ClosedLoopDynamics_velocity(this,x_k,v_desired)
            
            % Extract the current state properties
            p_k     = x_k(1:3);
            v_k     = x_k(4:6);
            R_k     = reshape(x_k(7:15),3,3); % The rotation matrix components
            eta_k   = OMAS_geometry.rotationMatrixToEulers(R_k);
            omega_k = x_k(16:18);
            
            % The state reference
            y_k = [p_k;v_k;eta_k;omega_k];
            
            % Extract the parameters from the state
            psi_desired = pi/2;
            y_desired = [zeros(3,1);v_desired;[0;0;psi_desired];zeros(3,1)];
            
            % Extract DYNAMIC parameters
            g  = this.DYNAMICS.g;
            m  = this.DYNAMICS.m;
            
            % State matrix
            A = this.DYNAMICS.A;
            A(4:6,7:9) = g*[ sin(psi_desired), cos(psi_desired), 0;
                            -cos(psi_desired), sin(psi_desired), 0;
                                            0,                0, 0];
            % Input matrix
            B = this.DYNAMICS.B;
            
            % LQR control gain
            [K,~,~] = lqr(A,B,eye(12),eye(4));
            K(:,1:3) = 0;
            % Calculate the error
            e_k  = y_k - y_desired;
            du_k = -K*e_k;
            % Calculate the inputs
            u_k = [m*g;0;0;0] + du_k;
            
            % Provide the inputs to the open-loop dynamics
            [dxdt] = this.OpenLoopDynamics(x_k,u_k);
        end
        % Quadcopter Dynamics (Closed-loop)
        function [dxdt] = ClosedLoopDynamics_position(this,x_k,p_desired)
            
            % Extract the current state properties
            p_k     = x_k(1:3);
            v_k     = x_k(4:6);
            R_k     = reshape(x_k(7:15),3,3); % The rotation matrix components
            eta_k   = OMAS_geometry.rotationMatrixToEulers(R_k);
            omega_k = x_k(16:18);
            
            % The state reference
            y_k = [p_k;v_k;eta_k;omega_k];
            
            % Extract the parameters from the state
            psi_desired = 0;
            y_desired = [p_desired(1:3);zeros(3,1);[0;0;psi_desired];zeros(3,1)];
            
            % Extract DYNAMIC parameters
            e3 = this.DYNAMICS.e3;
            g  = this.DYNAMICS.g;
            m  = this.DYNAMICS.m;
            I  = this.DYNAMICS.I;
            
            % State matrix
            A = this.DYNAMICS.A;
            A(1:3,4:6) = eye(3);
            A(4:6,7:9) = g*[ sin(psi_desired), cos(psi_desired), 0;
                -cos(psi_desired), sin(psi_desired), 0;
                0,                0, 0];
            A(7:9,10:12)=eye(3);
            
            % Input matrix
            B = this.DYNAMICS.B;
            B(4:6,1) = e3/m;
            B(10:12,2:4) = inv(I);
            
            % LQR control gain
            [K,~,~] = lqr(A,B,eye(12),eye(4));
            % Calculate the error
            e_k  = y_k - y_desired;
            du_k = -K*e_k;
            % Calculate the inputs
            u_k = [m*g;0;0;0] + du_k;
            
            % Provide the inputs to the open-loop dynamics
            [dxdt] = this.OpenLoopDynamics(x_k,u_k);
        end
        % Quadcopter Dynamics (Open-loop)
        function [dxdt] = OpenLoopDynamics(this,x_k,u_k)
            
            % Extract the state parameters
            v_k = x_k(4:6);
            R_k = reshape(x_k(7:15),3,3);
            w_k = x_k(16:18);
            % Extract the inputs
            f   = u_k(1);
            tau = u_k(2:4);
            
            % Get the constants
            m  = this.DYNAMICS.m;
            g  = this.DYNAMICS.g;
            I  = this.DYNAMICS.I;
            e3 = this.DYNAMICS.e3;
            
            % nonlinear dynamics based on rotation dynamics
            p_dot   = v_k;
            v_dot   = f/m*R_k*e3-g*e3;
            R_dot    = R_k*skew(w_k);
            vecR_dot = reshape(R_dot,9,1);
            w_dot    = inv(I)*(tau-skew(w_k)*I*w_k);
            
            % THE STATE DIFFERENCE
            dxdt = [p_dot;v_dot;vecR_dot;w_dot];
        end
    end
    % Dynamics container
    methods
        % Get the (generic) quadcopter dynamic properties
        function [DYNAMICS] = CreateDYNAMICS(this)
            % Simple quadcopter, use the default properties
            DYNAMICS = this.CreateDYNAMICS_default();
        end
        % Create (quadcopter) default dynamic structure
        function [DYNAMICS] = CreateDYNAMICS_default(this)
            % Define an infinite "throw" range by default
            DYNAMICS = struct();
            % Define kinematic constraints
            DYNAMICS.maxLinearVelocity      = ones(3,1)*this.v_max;	% Limits on the agents linear velocity
            DYNAMICS.maxLinearAcceleration  = inf(3,1);            	% Limits on the agents linear acceleration
            DYNAMICS.maxAngularVelocity     = ones(3,1)*this.w_max;	% Limits on the agents angular velocity
            DYNAMICS.maxAngularAcceleration = inf(3,1);            	% Limits on the agents angular acceleration
            % Dynamical properties
            DYNAMICS.e3 = [0;0;1];                                   % Local z-axis
            DYNAMICS.g  = 300978377846937375/30680772461461504;      % Gravitational constant
            DYNAMICS.m  = 1;                                         % Quadcopter mass
            DYNAMICS.I  = 2*diag([4.856, 4.856, 8.801])*10^(-3);     % Quadcopter inertia
            % Control properties
            % Plant matrix
            DYNAMICS.A = zeros(12);
            DYNAMICS.A(4:6,7:9)   = eye(3)*DYNAMICS.g;
            DYNAMICS.A(1:3,4:6)   = eye(3);
            DYNAMICS.A(7:9,10:12) = eye(3);
            % Input matrix
            DYNAMICS.B = zeros(12,4);
            DYNAMICS.B(4:6,1)     = DYNAMICS.e3/DYNAMICS.m;
            DYNAMICS.B(10:12,2:4) = inv(DYNAMICS.I);
            % Observation matrix
            DYNAMICS.C = eye(12);
            % Feed-forward matrix
            DYNAMICS.D = eye(4);
        end
    end
end

