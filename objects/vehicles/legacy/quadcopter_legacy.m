%% LEGACY QUADCOPTER (from code provided by Dr. Shiyu Zhao)
% In this agent, the basic function remain unaltered from his original
% code and is simply mapped to the OMAS architecture.

% Author: James A. Douthwaite

classdef quadcopter_legacy < quadcopter
    properties
        % GLOBAL.position   - The position of the object in the global axes
        % GLOBAL.velocity   - The velocity of the object in the global axes
        % GLOBAL.quaternion - The quaternion representing the earth to rotated body
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function [this] = quadcopter_legacy(varargin)
            
            % Call the super class
            this@quadcopter(varargin);                                      % Create the super class 'agent'
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end    
        % Setup - X = [pn;vn;Rvec;omega]
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % This function is called in order to build a local state
            % vector intended for aerospace control simulation (NED).
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z dx dy dz R1-R9 p q r]

            % THIS MODEL OPERATES IN THE GLOBAL FRAME
            % Properties can be taken directly from the VIRTUAL properties
            p0 = this.GetGLOBAL('position');   % True global position   
            v0 = this.GetGLOBAL('velocity');   % True global position   

            % GET THE ROTATION MATRIX fixed local axes -> rotated global pose
            % R0eval = OMAS_geometry.quaternionToRotationMatrix(qn0); % For validation
            R0 = OMAS_geometry.eulersToRotationMatrix(localXYZrotations);
            
            % BUILD THE GLOBAL STATE VECTOR     
            %v0     = R0'*localXYZVelocity;             % The global velocity
            vecR0   = reshape(R0',9,1);                 % Defines local vector to global vector
            omega0  = zeros(3,1);
            x0      = [p0;v0;vecR0;omega0];                      
            % ASSIGN THE LOCAL FRD STATE
            this.SetGLOBAL('priorState',x0);
            this.localState = x0;
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
%             % DEFAULT BEHAVIOUR
%             desiredSpeed = 1;
%             desiredHeadingVector = [1;0;0];
%             desiredVelocity = desiredHeadingVector*desiredSpeed;
%             % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
%             % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
%             [obj,~,~] = obj.getAgentUpdate(ENV,varargin{1});
            
            % Compute the control input
            [this] = this.Controller(ENV,this.localState);
            
            %fprintf('Object ID %s\n',num2str(obj.objectID));

            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            this.DATA.inputNames = {'p (rad/s)','q (rad/s)','r (rad/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(16:18);                        % Record the control inputs
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % Controller loop ((for quadcopter with a global state vector)
        function [this] = Controller(this,ENV,x_k)
            
            % //////// ///// UPDATE THE LOCAL STATE VECTOR ////////////////
            [x_k_plus] = this.UpdateLocalState(ENV,x_k);

            % GET THE LATEST ROTATION MATRIX
            R_k_plus = reshape(x_k_plus(7:15),3,3);                      % Is fixed axes to global
            q_k_plus = OMAS_geometry.rotationMatrixToQuaternion(R_k_plus');

            % Record new state
            this.localState = x_k_plus;
            
            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            this = this.GlobalUpdate_direct(...
                this.localState(1:3),...
                this.localState(4:6),...
                q_k_plus);     
        end
        % CONTINUOUS DYNAMIC LOOP
        function [out] = fcn_dynamicsAndControl(obj,in)
            pn=in(1:3);
            vn=in(4:6);
            R_BG=reshape(in(7:15),3,3);
            omega=in(16:18);
            
            % parameters
            g=9.8;
            I=eye(3);
            e3=I(:,3);
            m=2*0.5;%kg
            M=2*diag([4.856, 4.856, 8.801])*10^(-3);
            
            % control law
            %----------case 1: equilibrium input + noise
            % f=m*g+4*randn/2;
            % % tau=zeros(3,1)+[0.1 0 0]'; % will rotate and crash!
            % tau=zeros(3,1)+0.5*randn(3,1)/2; % will rotate and crash!
            
            %----------case 2: linear LQR based on the linearized model
            pn_desired=[5 5 0]';
            psi_desired=0;
            x_desireConst=[pn_desired;zeros(3,1);[0,0,psi_desired]';zeros(3,1)];
            
            A=zeros(12); % state matrix
            A(1:3,4:6)=eye(3);
            A(4:6,7:9)=g*[ sin(psi_desired), cos(psi_desired), 0;
                          -cos(psi_desired), sin(psi_desired), 0;
                                            0 0 0];
            A(7:9,10:12)=eye(3);
            
            B=zeros(12,4); % input matrix
            B(4:6,1)=e3/m;
            B(10:12,2:4)=inv(M);
            
            [K,~,~] = lqr(A,B,eye(12),eye(4)); % LQR control gain

            eta = obj.fcn_EulerFromRotation(R_BG);
            delta_x = [pn;vn;eta;omega]-x_desireConst;
            delta_u = -K*delta_x; % delta u=K*delta x
            u=[m*g;0;0;0]+delta_u;
            f=u(1);%+randn/10;
            tau=u(2:4);%+randn/10;
            
            % nonlinear dynamics based on rotation dynamics
            pn_dot=vn;
            vn_dot=f/m*R_BG*e3-g*e3;
            
            R_dot=R_BG*skew(omega);
            vecR_dot=reshape(R_dot,9,1);
            omega_dot=inv(M)*(tau-skew(omega)*M*omega);
            
            % THE STATE DIFFERENCE 
            out=[pn_dot;vn_dot;vecR_dot;omega_dot];
        end
    end
    % STATIC FUNCTIONS
    methods (Static)
        % IMPORT THE MODEL PROPERTIES
        function [A,B,K_lqr] = GetModelProperties(dt,g,e3,m,M,psi_desired)
            % PLANT MATRIX
            A=zeros(12);
            A(1:3,4:6) = eye(3);
            A(4:6,7:9) = g*[sin(psi_desired), cos(psi_desired), 0;
                           -cos(psi_desired), sin(psi_desired), 0;
                                           0,                0, 0];
            A(7:9,10:12) = eye(3);
            % INPUT MATRIX
            B=zeros(12,4); % input matrix
            B(4:6,1)=e3/m;
            B(10:12,2:4)=inv(M);
            
            % LQR CONTROLLER                       
            isdescrete = 1;
            if isdescrete
                C = eye(12);
                D = [zeros(8,4);eye(4)];
                SYS = ss(A,B,C,D);
                [SYSd] = c2d(SYS,dt,'zoh');
                
                Q = eye(12)*1E0;
                R = eye(4)*1E4;                
                [K_lqr] = dlqr(SYSd.A,SYSd.B,Q,R);
            else
                Q = eye(12);
                R = eye(4);
                [K_lqr,~,~] = lqr(A,B,Q,R); % LQR control gain
            end
        end
        % ROTATION MATRIX -> EULERS
        function [rho] = fcn_EulerFromRotation(R)
            % R is the rotation from body to world frame
            
            % assume -pi/2<tht<pi/2
            
            %R31=-sin(tht)
            tht=-asin(R(3,1));
            
            %R32=sin(phi)*cos(tht), R33=cos(phi)*cos(tht)
            phi=atan2(R(3,2),R(3,3));
            
            %R21=cos(tht)*sin(psi), R11=cos(tht)*cos(psi)
            psi=atan2(R(2,1),R(1,1));
            
            rho=[phi,tht,psi]';
        end
        % EULERS -> ROTATION MATRIX
        function [Rbl] = fcn_Euler2Rotation(phi,theta,psi)
            eulers = [phi, theta, psi];
            Rx=[1   0         0;
                0   cos(eulers(1))  sin(eulers(1));
                0  -sin(eulers(1))  cos(eulers(1))];
            Ry=[ cos(eulers(2))  0   -sin(eulers(2));
                0        1    0;
                sin(eulers(2))  0    cos(eulers(2))];
            Rz=[ cos(eulers(3))    sin(eulers(3))  0;
                -sin(eulers(3))    cos(eulers(3))  0;
                0           0         1];
            Rbl=Rx*Ry*Rz;                      %Rotation matrix: object -> camera frame   or NED to body
        end
        % SKEW MATRIX
        function [out] = fcn_skewSymmetric(x)
            out=[0   -x(3)  x(2);
                 x(3)    0 -x(1);
                -x(2) x(1)    0];
        end
    end
    % OMAS MAPPINGS
    methods 
        % GET THE STATE UPDATE (USING ODE45)
        function [X] = UpdateLocalState(obj,TIME,X0)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                X = X0;
                return
            else
                opts = odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2);
                [~,Xset] = ode45(@(t,X) obj.fcn_dynamicsAndControl(X),[0 TIME.dt],X0,opts);
                X = Xset(end,:)';
            end
        end    
    end
end




















