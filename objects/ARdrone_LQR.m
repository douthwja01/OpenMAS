%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_LQR < ARdrone
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    % LQR SPECIFIC PARAMETERS
    properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTOR METHOD
        function obj = ARdrone_LQR(varargin)
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@ARdrone(varargin);                                         % Create the super class 'agent'                  
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            
            % GET THE CLOSED LOOP LINEAR MODEL
            % This function overrides the linear plant model with the
            % closed loop plant
            [obj.DYNAMICS.K_lqr] = obj.getLQRFeedback(obj.DYNAMICS.SS);
            
            % DYNAMIC CONSTRAINTS
            obj.nominalSpeed = 1; 
            obj.maxSpeed = 1;
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);            
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
        function [obj] = main(obj,TIME,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % DEFAULT BEHAVIOUR 
            desiredHeadingVector = [1;0;0];
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,~] = obj.getAgentUpdate(varargin{1});
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
            end
            
            mapAngle = pi;
            XYZ2NED = [ 1             0              0;
                        0 cos(mapAngle) -sin(mapAngle);
                        0 sin(mapAngle)  cos(mapAngle)];
                    
            NEDVelocity = (XYZ2NED*desiredHeadingVector)*obj.nominalSpeed; 
            
            % CALCULATE THE NEW STATE REFERENCE
            X_desired = [zeros(3,1);zeros(2,1);obj.localState(6);NEDVelocity;zeros(3,1)];         % Create a state reference
            
            % /////////////////// AIRCRAFT DYNAMICS ///////////////////////
            [obj.localState] = obj.updateNonLinearPlant(TIME,obj.localState,X_desired);
%             [X] = obj.updateLinearPlant(TIME,obj.localState,X_desired)            
            
            % \\\\\\\\\\\\\\\\\\\\\ GLOBAL UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\
            obj = obj.updateGlobalProperties_ARdrone(dt,obj.localState);
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
%             obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
%             obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
%             obj.DATA.inputNames = {'p (rad/s)','q (rad/s)','r (rad/s)'};
            obj.DATA.inputNames = {'x','y','z','phi','theta','psi','dx','dy','dz','dphi','dtheta','dpsi'};
            obj.DATA.inputs(1:numel(obj.DATA.inputNames),TIME.currentStep) = obj.localState;          % Record the control inputs 
        end
    end
    
    % ///////////////////////// CONTROLLER ////////////////////////////////
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
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
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
            dU = - obj.DYNAMICS.K_lqr*Xerror;
            % CALL THE OPENLOOP DYNAMICS
            [dX] = obj.ARdrone_linear_openLoop(obj.DYNAMICS.SS,X,dU);
        end
    end
    % ///////////////////////// CONTROLLER ////////////////////////////////    
    methods (Static)    
        % GET LINEAR-QUADRATIC-REGULATOR(LQR) CONTROLLER
        function [K_lqr] = getLQRFeedback(SYS_openLoop)
            % This function designs the LQR controller used to provide
            % error feedback.
            
            % GET THE DESCRETE LQR CONTROLLER
            Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1.5E1;
            R = diag(ones(4,1))*1E-4;
            
            [K_lqr,~,~] = lqr(SYS_openLoop.A,SYS_openLoop.B,Q,R);  % Get the descrete LQR feedback  
            
            % IGNORE POSITION FEEDBACK
            K_lqr(:,1:3) = 0; % Omit position feedback
            % IGNORE XY POSITION FEEDBACK
%             K_lqr(:,1:2) = 0; % Omit XY position feedback
            % IGNORE HEADING FEEDBACK
%             K_lqr(:,6) = 0;  
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















