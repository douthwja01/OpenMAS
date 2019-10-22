%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_LQR < ARdrone
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    % LQR SPECIFIC PARAMETERS
    properties

    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function obj = ARdrone_LQR(varargin)
            
            % Call the super class
            obj@ARdrone(varargin);                                         % Create the super class 'agent'                  
             
            % DYNAMIC CONSTRAINTS
            obj.nominalSpeed = 2.0; 
            obj.maxSpeed = 3;
            
            % Append control parameter
            %obj.DYNAMICS.Q = diag([1 1 1 1 1 1 2 1 1 1 1 1])*1.5E1;
            obj.DYNAMICS.Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1.5E1;
            obj.DYNAMICS.R = diag(ones(4,1))*1E-4;     
            obj.DYNAMICS.N = zeros(size(obj.DYNAMICS.SS.B));        % Terminal input penalisation     
                        
            % Calculate the linear feedback from the LQR
            obj.DYNAMICS.K_lqr = obj.CreateController_LQR(...
                obj.DYNAMICS.SS,...
                obj.DYNAMICS.Q,...
                obj.DYNAMICS.R,...
                obj.DYNAMICS.N);
      
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
            desiredVelocity = desiredHeadingVector*obj.nominalSpeed;
            
            % Express velocity vector in NED coordinates
            if ENV.currentTime <= 1
                desiredVelocity = zeros(3,1);
            end
            
            % Convert from 'xyz' to 'ned'
            desiredVelocity = enu2ned(desiredVelocity);
            
            % LQR controller update
            obj = obj.controller(ENV,desiredVelocity);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                                   '$\dot{\phi}$ (rad/s)','$\dot{\theta}$ (rad/s)','$\dot{\psi}$ (rad/s)'};
            obj.DATA.inputs(1:numel(obj.DATA.inputNames),ENV.currentStep) = obj.localState(7:end);          % Record the control inputs 
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    % CONTROLLER METHODS
    methods
        % ARdrone controller function
        function [obj] = controller(obj,ENV,NEDVelocity)
            
            % CALCULATE THE NEW STATE REFERENCE
            if ~obj.IsIdle()
                X_desired = [zeros(5,1);obj.localState(6);NEDVelocity;zeros(3,1)];  
            else
                % Idle reference
                X_desired = [zeros(5,1);obj.localState(6);zeros(6,1)];
            end
            
            % /////////////////// AIRCRAFT DYNAMICS ///////////////////////
            useLinearModel = false;
            if ~useLinearModel
                [obj.localState] = obj.updateNonLinearPlant(ENV,obj.localState,X_desired);
            else
                [obj.localState] = obj.updateLinearPlant(ENV,obj.localState,X_desired);         
            end
            
            % \\\\\\\\\\\\\\\\\\\\\ GLOBAL UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\
            obj = obj.updateGlobalProperties_ARdrone(ENV.dt,obj.localState);
        end
    end
    methods (Static)
        % Create an LQR controller (feedback gain)
        function [K_lqr] = CreateController_LQR(SS,Q,R,N)
            % This function designs the LQR controller used to provide
            % error feedback.
            
            if nargin < 6
                N = zeros(size(SS.B));
            end
            
            % Get the LQR feedback             
            [K_lqr,~,~] = lqr(SS,Q,R,N);   
            
            % IGNORE POSITION FEEDBACK
            K_lqr(:,1:3) = 0;   % Omit position feedback
            % IGNORE XY POSITION FEEDBACK
%             K_lqr(:,1:2) = 0; % Omit XY position feedback
            % IGNORE HEADING FEEDBACK
%             K_lqr(:,6) = 0;  
        end 
    end
    % //////////////////////// MODELLING/DYNAMICS /////////////////////////
    methods
        % ODE45 - UPDATE NONLINEAR CLOSED LOOP DYNAMICS
        function [X] = updateNonLinearPlant(obj,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Input sanity check
            assert(isstruct(TIME),'Expecting a time structure.');
            
            X = X0; % Default to no change
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime ~= TIME.timeVector(end)
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                [~,Xset] = ode45(@(t,X) obj.ARdrone_nonLinear_closedLoop(X,X_desired),...
                    [0 TIME.dt],X0,odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)';
            end
        end
        % ODE45 - UPDATE LINEAR CLOSED LOOP DYNAMICS
        function [X] = updateLinearPlant(obj,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Input sanity check
            assert(isstruct(TIME),'Expecting a time structure.');
            
            X = X0; % Default to no change
            
            % Check for the last time-step
            if TIME.currentTime ~= TIME.timeVector(end) 
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP
                [~,Xset] = ode45(@(t,X) obj.ARdrone_linear_closedLoop(X,X_desired),...
                    [0 TIME.dt],X0,odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)';
            end
        end
    end
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
            [dX] = obj.ARdrone_linear_openLoop(X,dU);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















