%% SIMPLE FIXED WING AGENT (fixedWing.m) //////////////////////////////////
% This class is designed to contain the placeholder dynamics of a
% fixed-wing aircraft.

% Author: James A. Douthwaite

classdef fixedWing < agent
    properties

    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function [this] = fixedWing(varargin)
            % Call the super class
            this@agent(varargin);                                           % Create the super class 'agent'        
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup - 3D DUBLINS CAR STATE [x;y;z;v;psi;theta]
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % Infer the initial conditions from the scenario.
            p     = zeros(3,1);            % The position in the local frame
            speed = norm(localXYZVelocity);% The speed of the aircraft
            theta = localXYZrotations(2);  % The elevation of the aircraft
            psi   = 0;                     % The yaw of the aircraft
            % Define the initial state from the scenario
            this.localState = [p(1);p(2);p(3);speed;theta;psi]; 
            this = this.SetGLOBAL('priorState',this.localState); 
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            %[obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(varargin{1});                   % IDEAL INFORMATION UPDATE 
            
            u_k = [0.1;0.1;0.1]; % Acceleration, pitch rate, yaw rate
            x_k = this.localState;
            
            % UPDATE LOCAL STATE
            newState = this.UpdateLocalState(ENV,x_k,u_k);
            
            % UPDATE THE GLOBAL PROPERTIES OF THE AGENT
            [this] = this.GlobalUpdate(ENV.dt,newState);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            this.DATA.inputNames = {'$v (m/s)$','$\dot{\theta} (rad/s)$','$\dot{\psi} (rad/s)$'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = newState(4:6);         % Record the control inputs
        end
    end
    methods  
        % GET THE STATE UPDATE (USING ODE45)
        function [X] = UpdateLocalState(this,TIME,X0,u)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Is of the form: u = [v;a;psi_dot;theta_dot];
            X = X0;
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                return
            else
                [~,Xset] = ode45(@(t,X) this.dynamics_nonholonomicDublinsCar(X,u),...
                    [0 TIME.dt],X0,...
                    odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)'; % Pass the state at the end of the sample period
            end
        end
        % DEFINE THE GLOBAL UPDATE PROCEDURE FOR A 3D DUBLINS CAR MODEL
        function [this] = GlobalUpdate(this,dt,eulerState)
            
            % Retrieve current properties
            p_k	= this.GetGLOBAL('position');
            v_k	= this.GetGLOBAL('velocity');
            q_k	= this.GetGLOBAL('quaternion'); 
            x_k	= this.GetGLOBAL('priorState');
            
            % USE THE 'RATE' STATES DIRECTLY
            v_k_plus = [eulerState(4);0;0];                         % Assumed forward in the local axes
            omega_k_plus   = [0;eulerState(5);eulerState(6)];             % Neglect roll only
            localAxisRates = (omega_k_plus - [0;x_k(5:6)])/dt;
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % Previous properties
  
            % PREVIOUS ROTATION-MATRIX
            %R_k = quat2rotm(quaternion_k');
            R_k = OMAS_geometry.quaternionToRotationMatrix(q_k);
            % THE GLOBAL AXIS RATES       
            omega = R_k'*localAxisRates;
            % UPDATE THE QUATERNION POSE
            q_k_plus = OMAS_geometry.integrateQuaternion(q_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(q_k_plus');        
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            v_k_plus = R_k_plus'*v_k_plus;
            p_k_plus = p_k + dt*v_k;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.GlobalUpdate_direct(...
                p_k_plus,...    % Global position at k plius
                v_k_plus,...    % Global velocity at k plus
                q_k_plus);      % Quaternion at k plus
        end
    end
    % STATIC FUNCTIONS
    methods (Static)
        % THE NON-LINEAR HOLONOMIC DUBLINS-CAR MODEL
        function [dX] = dynamics_nonholonomicDublinsCar(X,u)
            % The state:  [dx;dy;dz;dv;dpsi;dtheta]
            % The inputs: [a,theta_dot,psi_dot]
%             dX = zeros(6,1);  % GLOBAL SPACE
%             dX(1) = X(4)*cos(X(6))*cos(X(5));
%             dX(2) = X(4)*cos(X(6))*sin(X(5));
%             dX(3) = X(4)*sin(X(6));
%             dX(4) = a;          % Acceleration
%             dX(5) = psi_dot;    % Yaw rate
%             dX(6) = theta_dot;  % Pitch rate
            dX = zeros(6,1);
            dX(1) = X(4);
            dX(2) = 0;
            dX(3) = 0;
            dX(4) = u(1);  % Acceleration
            dX(5) = u(2);  % Pitch rate
            dX(6) = u(3);  % Yaw rate
        end
    end
end