%% ARDRONE DYNAMIC MODEL (ARdrone.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class contains the dynamic properties of an ARdrone, along with all
% the necessary control parameters that can be derived from the model
% parameters.

% Author: James A. Douthwaite

classdef ARdrone < agent
%%% ARdrone CHILD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the ARdrone specific properties for importing a
    % small quadrotor into the environment.
    
    %% ////////////////// ADRONE SPECIFIC PARAMETERS //////////////////////
    properties
        % Performance Parameters
        battery_capacity  = 1000;     % Maximum battery capacity (mAh)
        battery_voltage   = 11.1;     % Maximum battery voltage (V)
        flight_time       = 720;      % Rated flight time (s)
        dt = -1;                      % The linearisation time-step (Assume continuous intially)
    end
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods 
        % CONSTRUCTOR
        function obj = ARdrone(varargin)
            % This function is to construct the ARdrone object using the
            % object defintions held in the 'agent' class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object

            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent(varargin);
            
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            
            % PARSE THE USER INPUTS FOR PARAMETER UPDATES
            [obj] = obj.configurationParser(obj,varargin);
            
            % GET THE ARDRONE DYNAMIC PROPERTIES
            DYNAMICS = obj.GetDyanamicProperties(obj.dt);
            obj.DYNAMICS = obj.configurationParser(DYNAMICS,varargin);
            
            % GET THE ARDRONE SENSOR PROPERTIES
            [AR_SENSORS] = obj.GetSensorProperties();
            sensorfields = fieldnames(AR_SENSORS);
            for i = 1:numel(sensorfields)
                obj.SENSORS.(sensorfields{i}) = AR_SENSORS.(sensorfields{i});
            end
            
            % SET THE RADIUS
            obj = obj.SetRadius(obj.DYNAMICS.length/2);                    % Body just used to determine the size parameter                       
        end
        % SETUP - LOCAL STATE [x;x_dot]'
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function generates the initial state vector for the
            % ARdrone model from the global parameters provided from the
            % scenario.
            
            eulerIndices = 4:6;
            velocityIndices = 7:9;
            
            % MAPPING FROM XYZ CONVENTION TO THE NED CONVENTION
            mapAngle = pi;
            XYZ2NED = [ 1             0              0;
                        0 cos(mapAngle) -sin(mapAngle);
                        0 sin(mapAngle)  cos(mapAngle)];
            
            % DEFINE THE INITIAL STATE
            X = zeros(12,1);
            X(eulerIndices) = XYZ2NED*localXYZrotations;
            X(velocityIndices) = XYZ2NED*localXYZVelocity;
            % APPEND TO VIRTUAL PROPERTIES
            obj.VIRTUAL.priorState = X;                                    % Record the previous state
            obj.localState = X;  
        end
        
        % THIS CLASS HAS NO 'main' CYCLE
        % The ARdrone class is designed to simply provide the dynamics of a
        % ARdrone quadrotor MAV.
        
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % GET THE STATE UPDATE (USING ODE45)
        function [X] = updateNonLinearPlant(obj,TIME,X0,U)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                X = X0;
                return
            else
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                opts = odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2);
                [~,Xset] = ode45(@(t,X) obj.ARdrone_nonLinear_openLoop(X,U),[0 TIME.dt],X0,opts);
                X = Xset(end,:)';
            end
        end 
        % GET THE STATE UPDATE (USING ODE45)
        function [X] = updateLinearPlant(obj,TIME,X0,U)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                X = X0;
                return
            else
                % INTEGRATE THE DYNAMICS OVER THE TIME STEP 
                opts = odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2);
                [~,Xset] = ode45(@(t,X) obj.ARdrone_linear_openLoop(obj.SYS,X,U),[0 TIME.dt],X0,opts);
                X = Xset(end,:)';
            end
        end
    end
    % //////////////////////// MODEL PROPERTIES ///////////////////////////
    methods (Static)
        % THE NONLINEAR AR-DRONE PLANT ODES
        function [dX] = ARdrone_nonLinear_openLoop(X,U)
            %#codegen
            
            %    This function was generated by the Symbolic Math Toolbox version 8.1.
            %    17-Jul-2018 15:13:41
            
            PHI_t = X(4,:);
            PHI_dot = X(10,:);
            PSI_dot = X(12,:);
            THETA_t = X(5,:);
            THETA_dot = X(11,:);
            omega_1 = U(1,:);
            omega_2 = U(2,:);
            omega_3 = U(3,:);
            omega_4 = U(4,:);
            x_dot = X(7,:);
            y_dot = X(8,:);
            z_dot = X(9,:);
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
            dX = [x_dot;y_dot;z_dot;PHI_dot;THETA_dot;PSI_dot;t5.*(-9.81e2./1.0e2)+THETA_dot.*t4.*y_dot+THETA_dot.*t2.*z_dot-PSI_dot.*t2.*t3.*y_dot+PSI_dot.*t3.*t4.*z_dot;-PHI_dot.*z_dot+t3.*t4.*(9.81e2./1.0e2)+PSI_dot.*t5.*z_dot-THETA_dot.*t4.*x_dot+PSI_dot.*t2.*t3.*x_dot;t6.*(-3.73475e-5)-t8.*3.73475e-5-t9.*3.73475e-5-t10.*3.73475e-5+PHI_dot.*y_dot+t2.*t3.*(9.81e2./1.0e2)-PSI_dot.*t5.*y_dot-THETA_dot.*t2.*x_dot-PSI_dot.*t3.*t4.*x_dot;t11+t12+THETA_dot.*omega_1.*3.231842180365297e-3-THETA_dot.*omega_2.*3.231842180365297e-3+THETA_dot.*omega_3.*3.231842180365297e-3-THETA_dot.*omega_4.*3.231842180365297e-3-t6.*t7.*8.516167950913242e-5-t7.*t8.*8.516167950913242e-5-THETA_dot.^2.*t4.*9.996789383561644e-1-PSI_dot.*THETA_dot.*t2-t3.*t4.*t13+PSI_dot.*THETA_dot.*t2.*t3.*9.996789383561644e-1;-t11+t12+PHI_dot.*PSI_dot-PHI_dot.*omega_1.*3.231842180365297e-3+PHI_dot.*omega_2.*3.231842180365297e-3-PHI_dot.*omega_3.*3.231842180365297e-3+PHI_dot.*omega_4.*3.231842180365297e-3+t6.*t7.*8.516167950913242e-5-t7.*t8.*8.516167950913242e-5-t5.*t13+PHI_dot.*THETA_dot.*t4.*9.996789383561644e-1-PHI_dot.*PSI_dot.*t2.*t3.*9.996789383561644e-1;t6.*4.175141847767905e-4-t8.*4.175141847767905e-4+t9.*4.175141847767905e-4-t10.*4.175141847767905e-4-PHI_dot.*THETA_dot.*1.000321164757521+PHI_dot.*THETA_dot.*t2.*1.000321164757521+PSI_dot.*THETA_dot.*t5.*1.000321164757521+PHI_dot.*PSI_dot.*t3.*t4.*1.000321164757521];
        end
        % THE LINEAR (SS) AR-DRONE PLANT
        function [dX] = ARdrone_linear_openLoop(SYS,X,U)
            % This function computes the state update using the derived
            % linear model of the ARdrone. The system is assumed linearised
            % around the hover condition, therefore the state-space model is
            % considered to be the change about this condition (i.e. the
            % inputs are delta_omega.^2
            
            % COMPUTE THE STATE UPDATE (i.e. Xdot = (A*x + B*du))
            PHI_t = X(4,:);
            PHI_dot = X(10,:);
            PSI_dot = X(12,:);
            THETA_t = X(5,:);
            THETA_dot = X(11,:);
            omega_1 = U(1,:);
            omega_2 = U(2,:);
            omega_3 = U(3,:);
            omega_4 = U(4,:);
            x_dot = X(7,:);
            y_dot = X(8,:);
            z_dot = X(9,:);
            t2 = omega_1.^2;
            t3 = sqrt(2.0);
            t4 = omega_2.^2;
            t5 = omega_3.^2;
            t6 = omega_4.^2;
            t7 = t3.*t5.*1.703233590182648e-4;
            t8 = t3.*t6.*1.703233590182648e-4;
            dX = [x_dot;y_dot;z_dot;PHI_dot;THETA_dot;PSI_dot;THETA_t.*(-9.81e2./1.0e2);PHI_t.*(9.81e2./1.0e2);t2.*(-7.4695e-5)-t4.*7.4695e-5-t5.*7.4695e-5-t6.*7.4695e-5;t7+t8-t2.*t3.*1.703233590182648e-4-t3.*t4.*1.703233590182648e-4;-t7+t8+t2.*t3.*1.703233590182648e-4-t3.*t4.*1.703233590182648e-4;t2.*8.35028369553581e-4-t4.*8.35028369553581e-4+t5.*8.35028369553581e-4-t6.*8.35028369553581e-4];

        end
        % GET THE NOMINAL THRUST INPUTS FROM CURRENT STATE
        function [Uss] = ARdrone_nominalInput(X)
            %#codegen

            %    This function was generated by the Symbolic Math Toolbox version 8.1.
            %    17-Jul-2018 15:14:14

            PHI_t = X(4,:);
            THETA_t = X(5,:);
            t2 = 3.875420890636849e8;
            t3 = cos(PHI_t);
            t4 = cos(THETA_t);
            t5 = t3.*t4;
            t6 = sqrt(t5);
            t7 = t2.*t6.*4.675627089909692e-7;
            Uss = [t7;t7;t7;t7];
        end
        % CONTROL MAPPING FOR ARDRONE
        function [setpoint] = ARdrone_controlMapping(throttle,roll,pitch,yaw)
            % This function maps the body axis torques onto the angular
            % rate set points of the motors (omega_(1-4))
            setpoint(1) = throttle - roll + pitch - yaw;
            setpoint(2) = throttle - roll - pitch + yaw;
            setpoint(3) = throttle + roll - pitch - yaw;
            setpoint(4) = throttle + roll + pitch + yaw; % This distributes the control inputs to the respective 
        end        
        % DYNAMIC PROPERTIES
        function [DYNAMICS] = GetDyanamicProperties(dt,config)
            % This function imports the complete set of dynamic properties
            % for the ARdrone 2.0 quadrotor UAV. 
            
            % DEFAULT CONFIGURATIONS
            if nargin == 1
                config = 'INDOOR';
            end
            % ALTER THE DYNAMICS BASED ON THE CONFIGURATION
            switch upper(config)
                case 'NONE'
                    DYNAMICS = struct('mass',0.366,'width',0.450,'length',0.290);
                case 'OUTDOOR'
                    % Arm length equiv : 0.3196
                    DYNAMICS = struct('mass',0.400,'width',0.452,'length',0.452);
                case 'INDOOR'
                    DYNAMICS = struct('mass',0.436,'width',0.515,'length',0.515);
                otherwise
                    error('Configuration not recognised');
            end 
            % INERTIAL MATRIX (x) configuration (FICTIONAL I MATRIX)
            DYNAMICS.I = 1.0E-06.*[2.8032E+04,0.0000E+00,0.0000E+00;
                                   0.0000E+00,2.8032E+04,0.0000E+00;
                                   0.0000E+00,0.0000E+00,2.8023E+04]; 
            % GET THE LINEAR STATE SPACE MODEL 
%             I_b11 = DYNAMICS.I(1,1);
%             I_b22 = DYNAMICS.I(2,2);
%             I_b33 = DYNAMICS.I(3,3);
%             m_b = 0.400;
%             kt = 1.4939E-05;
%             kh = 1.1700E-05;
%             l = 0.3196;
            g = 300978377846937375/30680772461461504;
            % PLANT MATRIX
            A = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
                 0, 0, 0, 0,-g, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, g, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
             
             % INPUT MATRIX
        B = [                                                     0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                                                                  0,                                                      0,                                                      0,                                                     0;
                             -1377879548585735/18446744073709551616,                 -1377879548585735/18446744073709551616,                 -1377879548585735/18446744073709551616,                -1377879548585735/18446744073709551616;
             -(5504628796600011325*2^(1/2))/32318695617139134431232, -(5504628796600011325*2^(1/2))/32318695617139134431232,  (5504628796600011325*2^(1/2))/32318695617139134431232, (5504628796600011325*2^(1/2))/32318695617139134431232;
              (5504628796600011325*2^(1/2))/32318695617139134431232, -(5504628796600011325*2^(1/2))/32318695617139134431232, -(5504628796600011325*2^(1/2))/32318695617139134431232, (5504628796600011325*2^(1/2))/32318695617139134431232;
                                   27848632988697/33350523172745984,                      -27848632988697/33350523172745984,                       27848632988697/33350523172745984,                     -27848632988697/33350523172745984];

 
            % DEFINE THE CONTINUOUS-TIME STATESPACE MODEL
            C = eye(12);                                                   % Build the observation matrix
            D = zeros(size(A,1),size(B,2));                                % Build the feed-forward matrix
            
            % CONVERT TO THE DESCRETE STATE-SPACE MODEL
            SYS = ss(A,B,C,D);       % Matlab continous statespace object
%             SYS = c2d(SYS,dt,'zoh'); % Convert the continuous statespace model to descrete. 
            
            % GET THE INVERTED DYNAMICS MATRIX (STEADY STATE CALCULATION)
            SSmatrix = [(eye(size(SYS.A,2))-SYS.A),-SYS.B;
                                             SYS.C, SYS.D];                % Defines the relationship between the dynamics and the stead state set-point
            DYNAMICS.SS = SYS;
            DYNAMICS.invSSmatrix = pinv(SSmatrix);
            
            % DEFINE THE NOMINAL CONTROL INPUT (HOVER INPUT)
            aNom = [-g;0;0;0];
            DYNAMICS.Uss = sqrt(inv(B(9:12,:))*aNom);
            
            % THE TRANSFER FUNCTION (omega -> thrust)
%             descrete_tf = tf(1,[0.108 1],dt); 
%             [A,B,C,D] = tf2ss(1,[0.108 1])
%             DYNAMICS.rotorSS = ss(A,B,C,D);
            
            % CONSTANT DYNAMIC PARAMETERS 
            DYNAMICS.omega_max = 750;                                      % Maximum angular rotor speed (rad/s)
            DYNAMICS.omega_min = 0;
        end
        % SENSOR PROPERTIES
        function [SENSORS]  = GetSensorProperties()
            % This function gets the ARdrone 2.0's sensor-specific data and
            % imports them to OMAS's obj.SENSOR structure.
            
            % Build the sensor description
            SENSORS = struct('ultrasound_freq',40E3,...     % Ultrasonic range-finder frequency (Hz)
                             'ultrasound_range',6,...       % Ultrasonic range-fidner range (m)
                             'camera_freq',30);             % Camera capture frequency (fps)
        end
    end
    % ///////////////////////// OMAS INTERFACES ///////////////////////////
    methods 
        % UNIQUE GLOBAL UPDATE FOR [x;x_dot]'
        function [obj] = updateGlobalProperties_ARdrone(obj,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y z phi theta psi dx dy dz dphi dtheta dpsi].
            
            positionIndices = 1:3;
            angleIndices = 4:6;
            
            % EQUIVALENT RATES
            velocity_k_plus = eulerState(7:9,1);
            localAxisRates  = eulerState(10:12,1);
%             velocity_k_plus = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt;
%             localAxisRates  = (eulerState(angleIndices) - obj.VIRTUAL.priorState(angleIndices))/dt; 
            
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isprop(obj,'targetWaypoint')
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj.VIRTUAL.idleStatus = logical(true);
                    velocity_k_plus = zeros(numel(positionIndices),1);     % Freeze the agent
                    localAxisRates  = zeros(numel(angleIndices),1);         % Freeze the agent
                end
            end
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = obj.VIRTUAL.quaternion;   % XYZ QUATERNION
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(quaternion_k');     % XYZ ROTATION
            
%             mapAngle = pi;
%             XYZ2NED = [ 1             0              0;
%                         0 cos(mapAngle) -sin(mapAngle);
%                         0 sin(mapAngle)  cos(mapAngle)];
            mapAngle = pi;        
            NED2XYZ = [ 1             0              0;
                        0 cos(mapAngle) -sin(mapAngle);
                        0 sin(mapAngle)  cos(mapAngle)];
                    
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*(NED2XYZ*localAxisRates);
            
            % UPDATE THE QUATERNION POSE
            quaternion_k_plus = OMAS_geometry.integrateQuaternion(quaternion_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(quaternion_k_plus');
            
            
            
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*(NED2XYZ*velocity_k_plus);
            
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end
    end
end