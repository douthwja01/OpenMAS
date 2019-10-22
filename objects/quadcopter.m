%% SIMPLE QUADCOPTER (adapted from the works of Dr. Shiyu Zhao) ///////////

% Author: James A. Douthwaite

classdef quadcopter < agent
    properties
        % globalPosition - The position of the object in the global axes
        % globalVelocity - The velocity of the object in the global axes
        % quaternion     - The quaternion representing the earth to rotated body
        
        % DYNAMICS - All the models parameters are held in the DYNAMICS
        %            container field.
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function obj = quadcopter(varargin)
            
            % Call the super class
            obj@agent(varargin);                                           % Create the super class 'agent'

            % Assign parameters
            [obj.DYNAMICS] = obj.importModelProperties();   
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end    
        % ///////////////////// SETUP FUNCTION ////////////////////////////
        % QUADCOPTER STATE VECTOR [pn;vn;etn;omega]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a quadrotor
            % object.

            % CALL THE STANDARD [x_t;x_dot} INITIALISER
            obj = obj.initialise_3DVelocities(localXYZVelocity,localXYZrotations);
        end
        % //////////////////// OBJEC MAIN CYCLE ///////////////////////////
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,~,~] = obj.GetAgentUpdate(ENV,varargin{1});                   % IDEAL INFORMATION UPDATE 
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Get desired heading
            desiredHeadingVector = obj.GetTargetHeading();
            desiredVelocity = desiredHeadingVector*obj.nominalSpeed;
                    
            % Express velocity vector in NED coordinates
            if ENV.currentTime <= 1
                desiredVelocity = zeros(3,1);
            end
            
            % Call the controller loop
            obj = obj.controller(ENV,desiredVelocity);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            obj.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                                   '$\dot{\phi}$ (rad/s)','$\dot{\theta}$ (rad/s)','$\dot{\psi}$ (rad/s)'};
            obj.DATA.inputs(1:numel(obj.DATA.inputNames),ENV.currentStep) = obj.localState(7:end);          % Record the control inputs 
        end
    end
    methods
        % The quadcopter velocity controller 
        function [obj] = controller(obj,ENV,v_k)
            
            % If idle, sit still
            if obj.IsIdle
                v_k = zeros(3,1);               
            end            
            % Convert local velocity to state reference
            X_desired = [zeros(5,1);obj.localState(6);v_k;zeros(3,1)];
            % UPDATE THE LOCAL STATE
            newState = obj.updateLocalState(ENV,obj.localState,X_desired);
            
            % UPDATE THE GLOBAL PROPERTIES OF THE AGENT
            [obj] = obj.updateGlobalProperties_3DVelocities(ENV.dt,newState);
        end        
        % CONTROL SELECTION
        function [obj,Xref] = GetControllerType(obj,type,dt)
            % This functional allows multiple types of LQR control to be
            % applied. 
            
            assert(ischar(type),'Please provide a control type.');
            
            % NO WAYPOINT IS OBSERVED, ASSUME IDLE
            if isempty(obj.targetWaypoint)
                Xref = obj.localState;
                Xref(4:5) = zeros(2,1);  % Level attitude
                Xref(7:12) = zeros(6,1); % No velocity
                return
            end
            
%             [R_phi] = eulerToRotationMatrix_roll(phi)
            
            switch upper(type)
                case 'POSITION'
                    % DESIGN SET POINT
                    desiredPosition = obj.targetWaypoint.position;
                    headingAngle = obj.localState(6);
                    Xref = [desiredPosition;zeros(2,1);headingAngle;zeros(3,1);zeros(3,1)]; % Velocity control
                case 'VELOCITY'
                    % MODIFY FEEDBACK
                    obj.DYNAMICS.K(1:4,1:3) = zeros(4,3);              % Omit the position feedback 
                    % DESIGN SET POINT
                    headingAngle = obj.localState(6);
                    headingVector = obj.GetTargetHeading();
                    desiredVelocity  = headingVector*obj.nominalSpeed; 
                    Xref = [zeros(3,1);zeros(2,1);headingAngle;desiredVelocity;zeros(3,1)]; % Velocity control
                case 'HEADING'
                    % MODIFY FEEDBACK
                    obj.DYNAMICS.K(1:4,1:3) = zeros(4,3);              % Omit the position feedback 
                    % GET VECTOR PROJECTIONS
                    headingVector = obj.GetTargetHeading();
                    [dlambda,dtheta] = obj.getVectorHeadingAngles([1;0;0],headingVector);
                    desiredVelocity = headingVector*obj.nominalSpeed;
                                  
                    headingRate = ((obj.localState(6) - dlambda) - obj.localState(6));
                    ascentRate = obj.nominalSpeed*sin(dtheta);
                    % DESIGN STATE REFERENCE
                    Xref = [zeros(3,1);...
                            [0;0;obj.localState(6)];...
                            [desiredVelocity(1);0;ascentRate];...
                            [0;0;headingRate]];   
                otherwise
                    error('Control type not recognised.');
            end            
        end
        % THE CLOSED-LOOP DYNAMICS
        function [dX] = closedLoopDynamics(obj,X,X_desired)
            % Get the next state from the previous state
            
            % THE DYNAMICAL CONSTANTSs
            g   = obj.DYNAMICS.g;       % Gravitational constant
            e3  = obj.DYNAMICS.e3;      % Body z-axes
            m   = obj.DYNAMICS.m;       % Body point mass (kg)  
            K   = obj.DYNAMICS.K;       % LQR feedback matrix
            eta = obj.localState(4:6);  % Current eulers

            % Represent gravity in the body axes
            R = OMAS_geometry.eulersToRotationMatrix(eta);
            T_b = R*(-g*m*e3);       	% Thrust in the body axes
            % COMPUTE THE ABSOLUTE INPUT
            X_error = X - X_desired;
            inputs = [T_b(3);0;0;0] - K*(X_error);
            % COMPUTE THE OPEN-LOOP DYNAMICS
            [dX] = obj.openLoopDynamics(X,inputs);
        end
        % THE OPEN-LOOP DYNAMICS
        function [dX] = openLoopDynamics(obj,X,inputs)
            % The inputs to the system are torques and forces directly.
            % inputs = [Tb;tau_r;tau_p;tau_y]
            
            Vb     = X(7:9,1);
            omegab = X(10:12,1);
            
            % DYNAMICS PROPERTIES
            g  = obj.DYNAMICS.g;
            e3 = obj.DYNAMICS.e3;
            m  = obj.DYNAMICS.m;
            M  = obj.DYNAMICS.M;
            R  = OMAS_geometry.eulersToRotationMatrix(X(4:6)); % Body to navigation frame
            
            % LINEAR & ANGULAR INPUTS
            Tb  = inputs(1);
            tau = inputs(2:4);
            
            % CALCULATE THE STATE DIFFERENCE
            dX = zeros(12,1);
            dX(1:3)   = Vb;
            dX(4:6)   = omegab;
            dX(7:9)   = Tb*(e3/m) + R'*(e3*m*g);
            dX(10:12) = M\(tau - OMAS_geometry.skew(omegab)*M*omegab);
        end
        % GET THE STATE UPDATE (USING ODE45)
        function [X]  = updateLocalState(obj,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                X = X0;
                return
            else
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                tspan = [0 TIME.dt];
                opts = odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2);
                [~,Xset] = ode45(@(t,X) obj.closedLoopDynamics(X,X_desired), tspan, X0, opts);
                X = Xset(end,:)';
            end
        end    
    end
    % STATIC FUNCTIONS
    methods (Static)
        % IMPORT THE MODEL PROPERTIES
        function [params] = importModelProperties()
            
            % BUILD THE DYNAMIC PROPERTIES OF THE AGENT
            g  = 9.81;
            I  = eye(3);
            e3 = I(:,3);
            m  = 2*0.5;    % kg
            M  = 2*diag([4.856, 4.856, 8.801])*10^(-3);
                        
            % PLANT MATRIX FOR STATE VECTOR 
            % X = [x y z phi theta psi dx dy dz dphi dtheta dpsi];
            A = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
                 0, 0, 0, 0,-g, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, g, 0, 0, 0, 0, 0, 0, 0, 0; % Gravity in XYZ
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            
            % INPUT MATRIX
            B = zeros(12,4);          % input matrix
            B(7:9,1) = e3/m;          % linear velocity terms
            B(10:12,2:4) = inv(M);    % angular velocity terms
            % OTHER STATE SPACE MATRICES
            C = eye(12);
            D = [zeros(8,4);eye(4)];
            % LQR CONTROLLER                       
            Q = eye(12);
            R = eye(4);
            % LQR CONTINUOUS FEEDBACK
            [K_lqr,~,~] = lqr(A,B,Q,R);        
            
            % RE-ORGANISE THE DYNAMIC PARAMETERS
            params = struct('g',g,...      % Gravitational constant
                                'e3',e3,...     % Body axis Z-vector
                                 'm',1,...      % Body mass
                                 'M',M,...      % Inertia tensor 
                                 'A',A,...      % Plant matrix
                                 'B',B,...      % Input matrix
                                 'C',C,...      % Observation matrix
                                 'D',D,...      % Feedforward matrix
                                 'Q',Q,...      % State penalisation matrix
                                 'R',R,...      % Input penalisation matrix
                                 'K',K_lqr);    % LQR feedback
        end
    end
end