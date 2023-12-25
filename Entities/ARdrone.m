
classdef ARdrone < quadcopter
    properties
        config = 'OUTDOOR'; % Indoor/outdoor dynamics
        % Performance Parameters
        battery_capacity  = 1000;	% Maximum battery capacity (mAh)
        battery_voltage   = 11.1;	% Maximum battery voltage (V)
        flight_time       = 720; 	% Rated flight time (s)
    end
    methods
        % Constructor
        function [this] = ARdrone(varargin)
            % Call the super class
            this@quadcopter(varargin);
            
            % Create the dynamics of a quadcopter
            this.SENSORS  = this.CreateSENSORS();     % Get the sensor parameters
            this.DYNAMICS = this.CreateDYNAMICS();
            this.GEOMETRY = this.CreateGEOMETRY();
            this.radius   = this.GEOMETRY.length/2;   % Use the "length" to define the radius
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a quadrotor
            % object.
            
            % For a state vector defined as x_k:
            % [x,y,z,phi,theta,psi,dx,dy,dz,wx,wy,wz]
            
            % THIS MODEL OPERATES IN THE GLOBAL FRAME
            %             p0 = this.GetGLOBAL('position');   % True global position
            %             v0 = this.GetGLOBAL('velocity');   % True global velocity
            %             % GET THE ROTATION MATRIX fixed local axes -> rotated global pose
            %             eta0 = OMAS_geometry.eulersToRotationMatrix(localXYZrotations);
            %             % Build an initial (global) state vector
            %             vecR0   = reshape(R0',9,1);        % Defines local vector to global vector
            %             omega0  = zeros(3,1);
            %             x0      = [p0;v0;vecR0;omega0];
            
            p0 = this.GetGLOBAL('position');   % True global position
            v0 = this.GetGLOBAL('velocity');   % True global velocity
            eta0   = localXYZrotations;
            omega0 = zeros(3,1);
            x0 = [p0;eta0;v0;omega0];
            
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
            desiredPosition = [1;1;0];
            
            % Call the controller loop
            this = this.Controller(ENV,desiredPosition);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            this.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                '$p$ (rad/s)','$q$ (rad/s)','$r$ (rad/s)'};
            % Record the control inputs
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(7:end);          
        end
    end
    %% ////////////////////// AUXILLARY METHODS ///////////////////////////
    methods
        % Get the state update (using ode45)
        function [x_k_plus] = UpdateLocalState(this,ENV,x_k,y_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Integrate across the time delta "dt"
            opts = odeset('RelTol',1e-2,'AbsTol',ENV.dt*1E-1);
            
            [~,Xset] = ode45(@(t,X) this.ClosedLoopDynamics_position(X,y_desired),[0 ENV.dt],x_k,opts);
            %            [~,Xset] = ode45(@(t,X) this.ClosedLoopDynamics_velocity(X,y_desired),[0 ENV.dt],x_k,opts);
            
            % Assign the integral state
            x_k_plus = Xset(end,:)';
        end
        % ARdrone Dynamics (Closed-loop)
        function [dxdt] = ClosedLoopDynamics_position(this,x_k,p_desired)
            
            % Extract the current state properties
            %             p_k   = x_k(1:3,1);
            %             eta_k = x_k(4:6,1);
            %             v_k   = x_k(7:9,1);
            %             w_k   = x_k(10:12,1);
            %             x_k = [p_desired(1:3);zeros(3,1);[0;0;psi_desired];zeros(3,1)];
            
            psi_desired = 0;
            % Define the desired state
            x_desired = [...
                p_desired(1:3);...
                [0;0;psi_desired];...
                zeros(3,1);...
                zeros(3,1)];
            
            % Extract DYNAMIC parameters
            e3 = this.DYNAMICS.e3;
            g  = this.DYNAMICS.g;
            m  = this.DYNAMICS.m;
            I  = this.DYNAMICS.I;
            
            % X = [dx dy dz wx wy wz ddx ddy ddz dwx dwy dwz]
            % State matrix
            A = this.DYNAMICS.A;
            A(1:6,7:12) = eye(6);
            A(7:9,4:6)  = g*[ sin(psi_desired), cos(psi_desired), 0;
                             -cos(psi_desired), sin(psi_desired), 0;
                                            0,                0, 0];
%             A_prev = zeros(12);
%             A_prev(1:3,4:6) = eye(3);
%             A_prev(4:6,7:9) = g*[ sin(psi_desired), cos(psi_desired), 0;
%                                  -cos(psi_desired), sin(psi_desired), 0;
%                                                  0,                0, 0];
%             A_prev(7:9,10:12)=eye(3);
            
            % Input matrix
            B = this.DYNAMICS.B;
            B(7:9,1) = e3/m;
            B(10:12,2:4) = inv(I);
            
            % LQR control gain
            [K,~,~] = lqr(A,B,eye(12),eye(4));
            % Calculate the error
            e_k  = x_k - x_desired;
            du_k = -K*e_k;
            % Calculate the inputs
            u_k = [m*g;0;0;0] + du_k;
            
            % Provide the inputs to the open-loop dynamics
            [dxdt] = this.OpenLoopDynamics(x_k,u_k);
        end
        % Quadcopter Dynamics (Open-loop)
        function [dxdt] = OpenLoopDynamics(this,x_k,u_k)
            % Designed for a state vector defined as:
            % x_k = [x,y,z,phi,theta,psi,dx,dy,dz,p,q,r]
            
            % Extract the state parameters
            eta_k = x_k(4:6,1);
            v_k   = x_k(7:9,1);
            w_k   = x_k(10:12,1);
            
            % Define the rotation matrix from the euler angles
%             R_k = reshape(x_k(7:15),3,3);
            R_k = OMAS_geometry.eulersToRotationMatrix(eta_k);

            % Extract the inputs
            f   = u_k(1);
            tau = u_k(2:4);
            
            % Get the constants
            m  = this.DYNAMICS.m;
            g  = this.DYNAMICS.g;
            I  = this.DYNAMICS.I;
            e3 = this.DYNAMICS.e3;
            
            % nonlinear dynamics based on rotation dynamics
            p_dot    = v_k;
            v_dot    = f/m*R_k*e3-g*e3;
            % Euler angles
            %R_dot    = R_k*skew(w_k);
            %eta_dot  = OMAS_geometry.rotationMatrixToEulers(R_dot);
            eta_dot = w_k;
            %vecR_dot = reshape(R_dot,9,1);
            
            w_dot    = inv(I)*(tau-skew(w_k)*I*w_k);
            
            % THE STATE DIFFERENCE
            dxdt = [p_dot;eta_dot;v_dot;w_dot];
        end
        % Global update from the new state
        function [this] = GlobalUpdate(this,x_k_plus)
            % This function updates the global structure from the new state
            % definition: [p_k;eta_k;v_k;omega_k]
            
            % Sanity Check
            assert(IsColumn(x_k_plus,12),"Expecting a 12DOF state vector.");
            
            % Define the rotation matrix from the euler angles
            R_k_plus = OMAS_geometry.eulersToRotationMatrix(x_k_plus(4:6,1));
            % Define the new quaternion from R
            q_k_plus = OMAS_geometry.rotationMatrixToQuaternion(R_k_plus');
            
            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            [this] = this.GlobalUpdate_direct(...
                x_k_plus(1:3),...   % The global position is within the state
                x_k_plus(7:9),...   % The global velocity is within the state
                q_k_plus);          % Append the global quaternion
            
            % Ensure the local state is re-assigned
            this.localState = x_k_plus;
        end
    end
    methods (Static)
        % Define the open-loop dynamics
        function [dxdt] = ARdrone_OpenLoopDynamics(x_k,u_k)
            %#codegen
            
            %    This function was generated by the Symbolic Math Toolbox version 8.1.
            %    17-Jul-2018 15:13:41
            
            PHI_t = x_k(4,:);
            PHI_dot = x_k(10,:);
            PSI_dot = x_k(12,:);
            THETA_t = x_k(5,:);
            THETA_dot = x_k(11,:);
            omega_1 = u_k(1,:);
            omega_2 = u_k(2,:);
            omega_3 = u_k(3,:);
            omega_4 = u_k(4,:);
            x_dot = x_k(7,:);
            y_dot = x_k(8,:);
            z_dot = x_k(9,:);
            t2 = cos(PHI_t);
            t3 = cos(THETA_t);
            t4 = sin(PHI_t);
            t5 = sin(THETA_t);
            t6 = omega_1.^2;
            t7 = sqrt(2.0);
            t8 = omega_2.^2;
            t9 = omega_3.^2;
            t10 = omega_4.^2;
            t11 = t7.*t9.*8.516167950913242e-5;
            t12 = t7.*t10.*8.516167950913242e-5;
            t13 = PSI_dot.^2;
            dxdt = [x_dot;y_dot;z_dot;PHI_dot;THETA_dot;PSI_dot;t5.*(-9.81e2./1.0e2)+THETA_dot.*t4.*y_dot+THETA_dot.*t2.*z_dot-PSI_dot.*t2.*t3.*y_dot+PSI_dot.*t3.*t4.*z_dot;-PHI_dot.*z_dot+t3.*t4.*(9.81e2./1.0e2)+PSI_dot.*t5.*z_dot-THETA_dot.*t4.*x_dot+PSI_dot.*t2.*t3.*x_dot;t6.*(-3.73475e-5)-t8.*3.73475e-5-t9.*3.73475e-5-t10.*3.73475e-5+PHI_dot.*y_dot+t2.*t3.*(9.81e2./1.0e2)-PSI_dot.*t5.*y_dot-THETA_dot.*t2.*x_dot-PSI_dot.*t3.*t4.*x_dot;t11+t12+THETA_dot.*omega_1.*3.231842180365297e-3-THETA_dot.*omega_2.*3.231842180365297e-3+THETA_dot.*omega_3.*3.231842180365297e-3-THETA_dot.*omega_4.*3.231842180365297e-3-t6.*t7.*8.516167950913242e-5-t7.*t8.*8.516167950913242e-5-THETA_dot.^2.*t4.*9.996789383561644e-1-PSI_dot.*THETA_dot.*t2-t3.*t4.*t13+PSI_dot.*THETA_dot.*t2.*t3.*9.996789383561644e-1;-t11+t12+PHI_dot.*PSI_dot-PHI_dot.*omega_1.*3.231842180365297e-3+PHI_dot.*omega_2.*3.231842180365297e-3-PHI_dot.*omega_3.*3.231842180365297e-3+PHI_dot.*omega_4.*3.231842180365297e-3+t6.*t7.*8.516167950913242e-5-t7.*t8.*8.516167950913242e-5-t5.*t13+PHI_dot.*THETA_dot.*t4.*9.996789383561644e-1-PHI_dot.*PSI_dot.*t2.*t3.*9.996789383561644e-1;t6.*4.175141847767905e-4-t8.*4.175141847767905e-4+t9.*4.175141847767905e-4-t10.*4.175141847767905e-4-PHI_dot.*THETA_dot.*1.000321164757521+PHI_dot.*THETA_dot.*t2.*1.000321164757521+PSI_dot.*THETA_dot.*t5.*1.000321164757521+PHI_dot.*PSI_dot.*t3.*t4.*1.000321164757521];
        end
    end
    % Properties
    methods
        % Get the (ARdrone) SENSOR structure
        function [SENSORS]  = CreateSENSORS(this)
            % This function gets the ARdrone 2.0's sensor-specific data and
            % imports them to OMAS's obj.SENSOR structure.
            
            % Get the default SENSOR structure
            SENSORS = this.SENSORS;
            SENSORS.ultrasound_freq = 40E3;     % Ultrasonic range-finder frequency (Hz)
            SENSORS.ultrasound_range = 6;       % Ultrasonic range-fidner range (m)
            SENSORS.camera_freq = 30;           % Camera capture frequency (fps)
        end
        % Get the (generic) quadcopter dynamic properties
        function [DYNAMICS] = CreateDYNAMICS(this)
            
            % Import the default properties of a simple quadcopter
            DYNAMICS = this.CreateDYNAMICS_default();
            
            % ALTER THE DYNAMICS BASED ON THE CONFIGURATION
            switch upper(this.config)
                case 'NONE'
                    DYNAMICS.m = 0.366;
                case 'OUTDOOR'
                    DYNAMICS.m = 0.400; % Arm length equiv : 0.3196
                case 'INDOOR'
                    DYNAMICS.m = 0.436;
                otherwise
                    error('Configuration not recognised');
            end
            % Assign the ARdrone properties
%             DYNAMICS.I = 1.0E-06.*[...
%                 2.8032E+04,0.0000E+00,0.0000E+00;
%                 0.0000E+00,2.8032E+04,0.0000E+00;
%                 0.0000E+00,0.0000E+00,2.8023E+04];
%             
            % Control properties
            % Plant matrix
            DYNAMICS.A = zeros(12); 
%             DYNAMICS.A(4:6,7:9)   = eye(3)*DYNAMICS.g;
%             DYNAMICS.A(1:3,4:6)   = eye(3);
%             DYNAMICS.A(7:9,10:12) = eye(3);
            % Input matrix
            DYNAMICS.B = zeros(12,4);
%             DYNAMICS.B(4:6,1)     = DYNAMICS.e3/DYNAMICS.m;
%             DYNAMICS.B(10:12,2:4) = inv(DYNAMICS.I);
            % Observation matrix
            DYNAMICS.C = eye(12);
            % Feed-forward matrix
            DYNAMICS.D = eye(4);    
        end
        % Get the (ARdrone) GEOMETRY structure
        function [GEOMETRY] = CreateGEOMETRY(this)
            
            % Parse the geometry associated with the object
            [GEOMETRY] = this.GetObjectGeometry(this);
            
            % ALTER THE DYNAMICS BASED ON THE CONFIGURATION
            switch upper(this.config)
                case 'NONE'
                    GEOMETRY.width  = 0.450;
                    GEOMETRY.length = 0.290;
                case 'OUTDOOR'
                    GEOMETRY.width  = 0.452;
                    GEOMETRY.length = 0.452; % Arm length equiv : 0.3196
                case 'INDOOR'
                    GEOMETRY.width  = 0.515;
                    GEOMETRY.length = 0.515;
                otherwise
                    error('Configuration not recognised');
            end
        end
    end
end