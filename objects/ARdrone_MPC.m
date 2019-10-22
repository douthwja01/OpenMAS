%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_MPC < ARdrone_LQR
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    % MPC SPECIFIC PARAMETERS
    properties
        % All ARdrone MPC parameters are within .DYNAMICS
        horizon = 10;
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [obj] = ARdrone_MPC(varargin)
            
            % Call the super class
            obj@ARdrone_LQR(varargin);                                     % Create the super class 'agent'                  
            
            % Build the MPC controller on top of the LQR
            MPC = obj.CreateController_MPC();
            fn  = fieldnames(MPC);
            for i = 1:length(fn)
                obj.DYNAMICS.(fn{i}) = MPC.(fn{i});
            end
            
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Main 
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project        
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            obj = obj.GetAgentUpdate(ENV,varargin{1});
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Get desired heading
            desiredHeadingVector = obj.GetTargetHeading();
            
            % Express velocity vector in NED coordinates
            NEDVelocity = enu2ned(desiredHeadingVector)*obj.nominalSpeed;
            
            % Ask the controller to achieve set point
            obj = obj.controller(ENV,NEDVelocity);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputNames = {'$\dot{x}$','$\dot{y}$','$\dot{z}$',...
                                   '$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$'};
            obj.DATA.inputs(1:numel(obj.DATA.inputNames),ENV.currentStep) = obj.localState(7:end);          % Record the control inputs 
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % ARdrone controller function
        function [obj] = controller(obj,ENV,NEDVelocity)
            
            % CALCULATE THE NEW STATE REFERENCE
            X_desired = [zeros(3,1);zeros(2,1);obj.localState(6);NEDVelocity;zeros(3,1)];  
            
            % /////////////////// AIRCRAFT DYNAMICS ///////////////////////
            useLinearModel = 1;
            if ~useLinearModel
                [obj.localState] = obj.updateNonLinearPlant(ENV,obj.localState,X_desired);
            else
                [obj.localState] = obj.updateLinearPlant(ENV,obj.localState,X_desired);         
            end
            
            % \\\\\\\\\\\\\\\\\\\\\ GLOBAL UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\
            obj = obj.updateGlobalProperties_ARdrone(ENV,obj.localState);
        end
        % GET LINEAR-QUADRATIC-REGULATOR(LQR) CONTROLLER
        function [MPC] = CreateController_MPC(obj)
            % This function calculates MPC controller constants and makes
            % them available for the controller updates.
            
            % Constants
            Ts = 0.1;                   % Descrete sampling period
            SS = obj.DYNAMICS.SS;       % For clarity
            Hp = obj.horizon;           % The prediction hsorizon
            Q  = obj.DYNAMICS.Q;        % The state feedback matrix
            R  = obj.DYNAMICS.R;        % The input feedback matrix
            Nq = zeros(size(obj.DYNAMICS.SS.A));        % The terminal feedback matrix
            
            % Convert to continous
            SS = c2d(SS,Ts,'zoh');      % Convert the SS model to descrete
            
            % GENERATE THE PREDICTION MATRICES
            n = size(SS.A,1);           % Rows in A (states)
            m = size(SS.B,2);           % Columns in B (inputs)
            Cd = [];                    % Stacked output matrix C
            Fd = zeros(n*Hp,n);         % F matrix, state*horizon by states
            Gd = zeros(n*Hp,m*(Hp-1));  % G matrix, state*horizon by (inputs*(horizon-1))
            
            % Define the prediction matrices
            for h = 1:Hp
                % THE PREDICTION MATRICES 
                Fd(n*(h-1)+(1:n),:) = SS.A^h; % The plant predictions, F
                for j = 1:h                   % Move horizontally
                    Gd(n*(h-1)+(1:n),m*(j-1)+(1:m)) = SS.A^(h-j)*SS.B; % The control predictions, G
                end
                % THE STACK OUTPUT CONVERSION, Cd
                Cd = blkdiag(Cd,SS.C); 
            end
            % DEFINE THE STACKED Q AND R MATRICES Qd, Rd
            Qd = kron(eye(Hp-1),Q);
            Qd = blkdiag(Qd,Nq);             % Building State prediction weighting matrices (for horizon)
            Rd = kron(eye(Hp),R);            % Building Input prediction weighting matrices (for horizon)
            % Fd - State predictions matrix
            % Gd - Control predictions matrix
            % Cd - Stack output matrix

            % DEFINE QUADRATIC COEFFICIENT
            theta = Cd*Gd;  
            % The (time-invariant) hessian matrix
            H = theta'*Qd*theta + Rd;        % The Hessian                   
            
            % The output selection matrix 
            Z = [eye(m),zeros(m,(m*(Hp-1)))];% The optimal output selection matrix [u* = Z.U*]

            % Steady state input calculation 
            % SS_state -> SS_input (SS_A*[x_ss;u_ss]' = [ 0, r ]
            T = [(eye(size(SS.A))-SS.A), -SS.B;
                                   SS.C,  SS.D];  
                               
            % Append control parameters
            MPC.Ts = Ts;           % The descrete sampling period
            MPC.Hp = Hp;           % Horizon
            MPC.Cd = Cd;           % The Stack observation matrix
            MPC.Qd = Qd;           % Stack plant weighting
            MPC.Rd = Rd;           % Stack control weighting
            MPC.Fd = Fd;           % State predictions matrix
            MPC.Gd = Gd;           % Control predictions matrix  
            MPC.Z = Z;             % Output selection matrix
            MPC.T = T;             % The matrix relationship between SS state & SS input
            MPC.Tinv = pinv(T);
            MPC.theta = theta;     % The jacobian
            MPC.H = H;             % The hessian
        end
    end
    %% /////////////////////// MODELLING/DYNAMICS /////////////////////////
    methods
        % ODE45 - UPDATE NONLINEAR CLOSED LOOP DYNAMICS
        function [X] = updateNonLinearPlant(obj,ENV,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            X = X0; % Default to no change
            
            % Check for the last time-step
            if ENV.currentTime ~= ENV.timeVector(end)               
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP
                [~,Xset] = ode45(@(t,X) obj.ARdrone_nonLinear_closedLoop(X,X_desired),...
                    [0 ENV.dt],X0,odeset('RelTol',1e-2,'AbsTol',ENV.dt*1E-2));
                X = Xset(end,:)';
            end
        end
        % ODE45 - UPDATE LINEAR CLOSED LOOP DYNAMICS
        function [X] = updateLinearPlant(obj,ENV,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            X = X0; % Default to no change
            
            % Check for the last time-step
            if ENV.currentTime ~= ENV.timeVector(end)  
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP
                [~,Xset] = ode45(@(t,X) obj.ARdrone_linear_closedLoop(X,X_desired),...
                    [0 ENV.dt],X0,odeset('RelTol',1e-2,'AbsTol',ENV.dt*1E-2));
                X = Xset(end,:)';
            end
        end
    end
    %% //////////////////// (MPC) CLOSED-LOOP SYSTEMS /////////////////////
    methods
        % CLOSED-LOOP NONLINEAR DYNAMICS
        function [dX] = ARdrone_nonLinear_closedLoop(obj,X,X_desired)
            % THE STATE ERROR
            Xerror = X - X_desired;
            % THE NOMINAL INPUT
            [Uss] = obj.ARdrone_nominalInput(X);
            % GET THE NOMINAL INPUT + LQR FEEDBACK
            U = Uss - obj.DYNAMICS.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = obj.ARdrone_nonLinear_openLoop(X,U);
        end
        % CLOSED-LOOP LINEAR DYNAMICS
        function [dX] = ARdrone_linear_closedLoop(obj,X,X_desired)
            % THE STATE ERROR
            Xerror = X - X_desired;
            % GET THE NOMINAL INPUT + LQR FEEDBACK
            dU = - obj.DYNAMICS.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = obj.ARdrone_linear_openLoop(obj.DYNAMICS.SS,X,dU);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















