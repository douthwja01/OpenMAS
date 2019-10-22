%% SIMPLE FIXED WING AGENT (fixedWing.m) //////////////////////////////////
% This class is designed to contain the placeholder dynamics of a
% fixed-wing aircraft.

% Author: James A. Douthwaite

classdef fixedWing < agent
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
        function obj = fixedWing(varargin)
            % Call the super class
            obj@agent(varargin);                                           % Create the super class 'agent'
                        
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup - 3D DUBLINS CAR STATE [x;y;z;v;psi;theta]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % Infer the initial conditions from the scenario.
            p     = zeros(3,1);            % The position in the local frame
            speed = norm(localXYZVelocity);% The speed of the aircraft
            theta = localXYZrotations(2);  % The elevation of the aircraft
            psi   = 0;                     % The yaw of the aircraft
            % Define the initial state from the scenario
            obj.localState = [p(1);p(2);p(3);speed;theta;psi]; 
            obj = obj.SetVIRTUALparameter('priorState',obj.localState); 
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            %[obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(varargin{1});                   % IDEAL INFORMATION UPDATE 
            
            u = [0.1;0.1;0.1]; % Acceleration, pitch rate, yaw rate
            
            % UPDATE LOCAL STATE
            newState = obj.updateLocalState(ENV,obj.localState,u);
            
            % UPDATE THE GLOBAL PROPERTIES OF THE AGENT
            [obj] = obj.updateGlobalProperties(dt,newState);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            obj.DATA.inputNames = {'$v (m/s)$','$\dot{\theta} (rad/s)$','$\dot{\psi} (rad/s)$'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = newState(4:6);         % Record the control inputs
        end
    end
    methods        
        % DEFINE THE GLOBAL UPDATE PROCEDURE FOR A 3D DUBLINS CAR MODEL
        function [obj] = updateGlobalProperties(obj,dt,DCMState)
            % USE THE 'RATE' STATES DIRECTLY
            velocity_k_plus = [DCMState(4);0;0];                           % Assumed forward in the local axes
            eulers_k_plus   = [0;DCMState(5);DCMState(6)];                 % Neglect roll only
            localAxisRates  = (eulers_k_plus - [0;obj.VIRTUAL.priorState(5:6)])/dt;
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = obj.VIRTUAL.quaternion;   
            % PREVIOUS ROTATION-MATRIX
            %R_k = quat2rotm(quaternion_k');
            R_k = OMAS_geometry.quaternionToRotationMatrix(quaternion_k);
            % THE GLOBAL AXIS RATES       
            omega = R_k'*localAxisRates;
            % UPDATE THE QUATERNION POSE
            quaternion_k_plus = OMAS_geometry.integrateQuaternion(quaternion_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(quaternion_k_plus');        
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      DCMState);                % The new state for reference
        end
        % GET THE STATE UPDATE (USING ODE45)
        function [X]   = updateLocalState(obj,TIME,X0,u)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % Is of the form: u = [v;a;psi_dot;theta_dot];
            X = X0;
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                return
            else
                [~,Xset] = ode45(@(t,X) obj.dynamics_nonholonomicDublinsCar(X,u),...
                    [0 TIME.dt],X0,...
                    odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                X = Xset(end,:)'; % Pass the state at the end of the sample period
            end
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