%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_LQR < ARdrone
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    % LQR SPECIFIC PARAMETERS
    properties
        % GET THE DESCRETE LQR CONTROLLER
        Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1.5E1;
        R = diag(ones(4,1))*1E-4;
        K_lqr;
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % CONSTRUCTOR
        function obj = ARdrone_LQR(varargin)
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@ARdrone(varargin);                                         % Create the super class 'agent'                  
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            
            % GET THE CLOSED LOOP LINEAR MODEL
            % Calculate the linear feedback from the LQR
            obj.K_lqr = obj.GetLQRController(obj.DYNAMICS.SS);
            
            % DYNAMIC CONSTRAINTS
            obj.nominalSpeed = 1; 
            obj.maxSpeed = 1;
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);            
        end
        % MAIN 
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
            NEDVelocity = obj.GetNEDFromENU(desiredHeadingVector)*obj.nominalSpeed;
            
            % LQR controller update
            obj = obj.controller(ENV,NEDVelocity);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputNames = {'x','y','z','phi','theta','psi','dx','dy','dz','dphi','dtheta','dpsi'};
            obj.DATA.inputs(1:numel(obj.DATA.inputNames),ENV.currentStep) = obj.localState;          % Record the control inputs 
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    % CONTROLLER METHODS
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
            obj = obj.updateGlobalProperties_ARdrone(ENV.dt,obj.localState);
        end
        % GET LINEAR-QUADRATIC-REGULATOR(LQR) CONTROLLER
        function [K_lqr] = GetLQRController(obj,SS)
            % This function designs the LQR controller used to provide
            % error feedback.
            
            % GET THE DESCRETE LQR CONTROLLER
%             Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1.5E1;
%             R = diag(ones(4,1))*1E-4;
            
            [K_lqr,~,~] = lqr(SS.A,SS.B,obj.Q,obj.R);  % Get the descrete LQR feedback  
            
            % IGNORE POSITION FEEDBACK
            K_lqr(:,1:3) = 0; % Omit position feedback
            % IGNORE XY POSITION FEEDBACK
%             K_lqr(:,1:2) = 0; % Omit XY position feedback
            % IGNORE HEADING FEEDBACK
%             K_lqr(:,6) = 0;  
        end
    end
    
    methods
        % 0DE45 - UPDATE NONLINEAR CLOSED LOOP DYNAMICS
        function [X] = updateNonLinearPlant(obj,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                X = X0;
                return
            else
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                opts = odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2);
                [~,Xset] = ode45(@(t,X) obj.ARdrone_nonLinear_closedLoop(X,X_desired),[0 TIME.dt],X0,opts);
                X = Xset(end,:)';
            end
        end
        % 0DE45 - UPDATE LINEAR CLOSED LOOP DYNAMICS
        function [X] = updateLinearPlant(obj,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            assert(isstruct(TIME),'Expecting a time structure.');
            
            % DETERMINE THE INTEGRATION PERIOD
            %if TIME.currentTime == TIME.timeVector(end)
            if obj.isIdle()
                X = X0;
                return
            end
            
            % INTEGRATE THE DYNAMICS OVER THE TIME STEP
            [~,Xset] = ode45(@(t,X) obj.ARdrone_linear_closedLoop(X,X_desired),...
                [0 TIME.dt],...
                X0,...
                odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
            X = Xset(end,:)';
        end
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
            dU = - obj.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = obj.ARdrone_linear_openLoop(X,dU);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















