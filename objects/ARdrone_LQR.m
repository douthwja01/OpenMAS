%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_LQR < ARdrone_prev
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    % LQR SPECIFIC PARAMETERS
    properties

    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function this = ARdrone_LQR(varargin)
            
            % Call the super class
            this@ARdrone_prev(varargin);                                         % Create the super class 'agent'                  
             
            % DYNAMIC CONSTRAINTS
            this.v_nominal = 2.0; 
            this.v_max = 3;
            
            % Append control parameter
            %obj.DYNAMICS.Q = diag([1 1 1 1 1 1 2 1 1 1 1 1])*1.5E1;
            this.DYNAMICS.Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1.5E1;
            this.DYNAMICS.R = diag(ones(4,1))*1E-4;     
            this.DYNAMICS.N = zeros(size(this.DYNAMICS.SS.B));        % Terminal input penalisation     
                        
            % Calculate the linear feedback from the LQR
            this.DYNAMICS.K_lqr = this.CreateController_LQR(...
                this.DYNAMICS.SS,...
                this.DYNAMICS.Q,...
                this.DYNAMICS.R,...
                this.DYNAMICS.N);
      
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides 
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        % - The same as any 2D/3D object
        % Main 
        function this = main(this,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project

            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            this = this.GetAgentUpdate(ENV,varargin{1});
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Get desired heading
            desiredHeadingVector = this.GetTargetHeading();
            desiredVelocity = desiredHeadingVector*this.v_nominal;
            
            % Express velocity vector in NED coordinates
            if ENV.currentTime <= 1
                desiredVelocity = zeros(3,1);
            end
            
            p_wp = this.targetWaypoint.position(:,this.targetWaypoint.sampleNum);
            
            % Convert from 'xyz' to 'ned'
%             desiredVelocity = enu2ned(desiredVelocity);
%             desiredNEDstate = [p_wp;zeros(3,1);zeros(6,1)];
            
            % LQR controller update
            this = this.Controller_position(ENV,enu2ned(p_wp));
            
%             this = this.Controller_velocity(ENV,enu2ned(desiredVelocity));
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            this.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                                   '$\dot{\phi}$ (rad/s)','$\dot{\theta}$ (rad/s)','$\dot{\psi}$ (rad/s)'};
            this.DATA.inputs(1:numel(this.DATA.inputNames),ENV.currentStep) = this.localState(7:end);          % Record the control inputs 
        end
    end
    
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    % CONTROLLER METHODS
    methods
        % ARdrone controller (local position)
        function [this] = Controller_position(this,ENV,NEDposition)
            
            % CALCULATE THE NEW STATE REFERENCE
            if ~this.IsIdle()
                X_desired = [NEDposition;zeros(2,1);-this.localState(6);zeros(6,1)];  
            else
                % Idle reference
                X_desired = [zeros(5,1);this.localState(6);zeros(6,1)];
            end
            
            % /////////////////// AIRCRAFT DYNAMICS ///////////////////////
            useLinearModel = false;
            if ~useLinearModel
                [this.localState] = this.UpdateNonLinearPlant(ENV,this.localState,X_desired);
            else
                [this.localState] = this.UpdateLinearPlant(ENV,this.localState,X_desired);         
            end
            
            % \\\\\\\\\\\\\\\\\\\\\ GLOBAL UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\
            this = this.GlobalUpdate_ARdrone(ENV.dt,this.localState);
        end
        % ARdrone controller (local velocity)
        function [this] = Controller_velocity(this,ENV,NEDVelocity)
            
            % CALCULATE THE NEW STATE REFERENCE
            if ~this.IsIdle()
                X_desired = [zeros(5,1);this.localState(6);NEDVelocity;zeros(3,1)];  
            else
                % Idle reference
                X_desired = [zeros(5,1);this.localState(6);zeros(6,1)];
            end
            
            % /////////////////// AIRCRAFT DYNAMICS ///////////////////////
            useLinearModel = false;
            if ~useLinearModel
                [this.localState] = this.UpdateNonLinearPlant(ENV,this.localState,X_desired);
            else
                [this.localState] = this.UpdateLinearPlant(ENV,this.localState,X_desired);         
            end
            
            % \\\\\\\\\\\\\\\\\\\\\ GLOBAL UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\
            this = this.GlobalUpdate_ARdrone(ENV.dt,this.localState);
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
        function [X] = UpdateNonLinearPlant(this,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Input sanity check
            assert(isstruct(TIME),'Expecting a time structure.');
            
            X = X0; % Default to no change
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime ~= TIME.timeVector(end)
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                [~,Xset] = ode45(@(t,X) this.ARdrone_nonLinear_closedLoop(X,X_desired),...
                    [0 TIME.dt],X0,odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)';
            end
        end
        % ODE45 - UPDATE LINEAR CLOSED LOOP DYNAMICS
        function [X] = updateLinearPlant(this,TIME,X0,X_desired)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Input sanity check
            assert(isstruct(TIME),'Expecting a time structure.');
            
            X = X0; % Default to no change
            
            % Check for the last time-step
            if TIME.currentTime ~= TIME.timeVector(end) 
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP
                [~,Xset] = ode45(@(t,X) this.ARdrone_linear_closedLoop(X,X_desired),...
                    [0 TIME.dt],X0,odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)';
            end
        end
    end
    methods
        % CLOSED-LOOP NONLINEAR DYNAMICS
        function [dX] = ARdrone_nonLinear_closedLoop(this,X,X_desired)
            % THE STATE ERROR
            Xerror = X - X_desired;
            % THE NOMINAL INPUT
            [Uss] = this.ARdrone_nominalInput(X);
            % GET THE NOMINAL INPUT + LQR FEEDBACK
            U = Uss - this.DYNAMICS.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = this.ARdrone_nonLinear_openLoop(X,U);
        end
        % CLOSED-LOOP LINEAR DYNAMICS
        function [dX] = ARdrone_linear_closedLoop(this,X,X_desired)
            % THE STATE ERROR
            Xerror = X - X_desired;
            % GET THE NOMINAL INPUT + LQR FEEDBACK
            dU = - this.DYNAMICS.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = this.ARdrone_linear_openLoop(X,dU);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















