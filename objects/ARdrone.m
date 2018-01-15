%% ARDRONE DYNAMIC MODEL (ARdrone.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class contains the dynamic properties of an ARdrone, along with all
% the necessary control parameters that can be derived from the model
% parameters.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef ARdrone < agent
%%% ARdrone CHILD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the ARdrone specific properties for importing a
    % small quadrotor into the environment.
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties
        % Phyiscal Parameters are initialised in the constructor
        
        % Performance Parameters
        omega_max         = 750;      % Maximum angular rotor speed (rad/s)
        battery_capacity  = 1000;     % Maximum battery capacity (mAh)
        battery_voltage   = 11.1;     % Maximum battery voltage (V)
        flight_time       = 720;      % Rated flight time (s)
        ultrasound_freq   = 40E3;     % Ultrasonic range-finder frequency (Hz)
        ultrasound_range  = 6;        % Ultrasonic range-fidner range (m)
        camera_freq       = 30;       % Camera capture frequency (fps)
        % STATESPACE
        SYS;                          % The descrete state-space model
        invertedSYS;                  % The steady state-inverted plant parameters
        dt = -1;
        % ROTOR TRANSFER FUNCTION
        rotorTF;
        rotorSS;
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = ARdrone(varargin)
            % This function is to construct the ARdrone object using the
            % object defintions held in the 'agent' class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object

            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent(varargin);
            % CHECK FOR USER OVERRIDES BEFORE CONTINUING
            [obj] = obj.configurationParser(varargin);
            
            % APPEND THE STATESPACE PARAMETERS
            [obj.SYS,obj.invertedSYS] = obj.getLocalLinearProperties(obj.dt);
            
            % GET THE ROTOR TRANSFER FUNCTION
            [obj.rotorTF,obj.rotorSS] = obj.getRotorTransferFunction(obj.dt);
            
            % PHYISCAL PARAMETERS
            [BODY] = obj.getBodyParameters('INDOOR');
            obj.VIRTUAL.size = BODY.length;                                % Body just used to determine the size parameter        
        end
        
        % UPDATE ARDRONE STATE VECTOR [x;y;z;psi;the;phi;u;v;w;p;q;r]
        function [newState] = stateDynamics(obj,X0,dt,inputs)
            % This function evaluates the ARdrones new state given the
            % control inputs and current state.
            % INPUTS:
            % X0          - The previous object state vector
            % dt          - The unit time step
            % omega_(1-4) - The rotor angular velocities (rad/s)
            
            % APPLY ROTOR SPEED SATURATIONS
           inputs = obj.boundValue(inputs,0,obj.omega_max);
            
            % ROTOR TRANSFER FUNCTION
%             inputs = obj.rotorSS.A*inputs;
            
            % GET THE STATE UPDATE FROM THE DYNAMIC EXPRESSIONS
            Xdot = obj.ARdrone_nonlinear_local(X0,inputs(1),inputs(2),inputs(3),inputs(4));
            % GET THE NEW UPDATED STATE
            newState = X0 + dt*Xdot;                           % Access the current state and append difference
        end
    end
    methods (Static)
        %% CONTROL FUNCTIONS
        % CONTROL MAPPING FOR ARDRONE
        function [V_setpoint] = controlMapping_Xtype(throttle,roll,pitch,yaw)
            % This function maps the body axis torques onto the angular
            % rate set points of the motors (omega_(1-4))
            V_setpoint(1) = throttle - roll + pitch - yaw;
            V_setpoint(2) = throttle - roll - pitch + yaw;
            V_setpoint(3) = throttle + roll - pitch - yaw;
            V_setpoint(4) = throttle + roll + pitch + yaw; % This distributes the control inputs to the respective 
        end
        
        %% PLANT FUNCTIONS
        % GET ROTOR TRANSFER FUNCTION
        function [descrete_tf,SYS] = getRotorTransferFunction(dt)
           % The calculated transfer function was found to be:
           % omega_out = 1/[0.108s + 1)*omega_in
           
           % Define the descrete transfer function
           descrete_tf = tf(1,[0.108 1],dt); 
           [A,B,C,D] = tf2ss(1,[0.108 1])
           SYS = ss(A,B,C,D);
        end
        % INVERTED KINEMATIC MATRIX
        function [inputs] = inverseKinematics(desiredAccelerations)
            % This function calculates the inputs required to generate the
            % desired body accelerations ([x,phi,theta,psi])
            
            % CONSTANT INVERTED PLANT MATRIX
            inverseKinematics= 1.0e+04 *[-1.6735   -7.3439    7.3439   -2.1368;...
                                         -1.6735   -7.3439   -7.3439    2.1368;...
                                         -1.6735    7.3439   -7.3439   -2.1368;...
                                         -1.6735    7.3439    7.3439    2.1368];
            % SOLVE FOR THE INPUTS REQUIRED FOR THE BODY ACCELERATIONS
            inputs = inverseKinematics*desiredAccelerations;
            inputs = sqrt(inputs); % The reverse kinematics are in terms of omega^2
        end
        % THE LINEARISED PLANT STATESPACE MODEL (LOCAL)
        function [SYS,invertedSYS] = getLocalLinearProperties(dt)
            % This function generates the state-space model representing
            % the ARdrone dynamics linearised around the hover condition.
            fprintf('Getting the statespace model\n');
            % PLANT MATRIX
            A = [0, 0, 0,                                    0,                                     0, 0, 1, 0, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 1, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 1, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 1, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 1, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 0, 1;
                 0, 0, 0,                                    0, -300978377846937375/30680772461461504, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 300978377846937375/30680772461461504,                                     0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0,                                    0,                                     0, 0, 0, 0, 0, 0, 0, 0];
             % INPUT MATRIX
             B = [                                               0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                                                                 0,                                               0,                                               0,                                              0;
                          -34446988714643375/502673776008585281536,        -34446988714643375/502673776008585281536,        -34446988714643375/502673776008585281536,       -34446988714643375/502673776008585281536;
                   -125592194882400484375/517099129874226150899712, -125592194882400484375/517099129874226150899712,  125592194882400484375/517099129874226150899712, 125592194882400484375/517099129874226150899712;
                    125592194882400484375/517099129874226150899712, -125592194882400484375/517099129874226150899712, -125592194882400484375/517099129874226150899712, 125592194882400484375/517099129874226150899712;
                                 -27848632988697/33350523172745984,                27848632988697/33350523172745984,               -27848632988697/33350523172745984,               27848632988697/33350523172745984];
 
 
            
            % DEFINE THE CONTINUOUS-TIME STATESPACE MODEL
%             A = zeros(12);
%             A(1:6,7:12) = eye(6);
%             A(8,4) = 9.81;
%             A(7,5) = -9.81;                                                % Build the plant matrix
%             B = zeros(12,4);
%             B(9:12,1:4) = [-0.0685,-0.0685,-0.0685,-0.0685;
%                            -0.2429,-0.2429, 0.2429, 0.2429;
%                             0.2429,-0.2429,-0.2429, 0.2429;
%                            -0.8350, 0.8350,-0.8350, 0.8350]*1E-3;          % Build the input matrix
            C = eye(12);                                                   % Build the observation matrix
            D = zeros(size(A,1),size(B,2));                                % Build the feed-forward matrix
            
            % CONVERT TO THE DESCRETE STATE-SPACE MODEL
            SYS = ss(A,B,C,D);                                             % Matlab continous statespace object
            SYS = c2d(SYS,dt);                                             % Convert to discrete      
%             SYS = ss(A,B,C,D);                                           % Matlab statespace object
            % GET THE INVERTED DYNAMICS MATRIX (STEADY STATE CALCULATION)
            SSmatrix = [(eye(size(SYS.A,2))-SYS.A),-SYS.B;
                                             SYS.C, SYS.D];                % Defines the relationship between the dynamics and the stead state set-point
            invertedSYS = pinv(SSmatrix);
        end
        % THE NONLINEAR AR-DRONE PLANT ODES (LOCAL)
        function [dX] = ARdrone_nonlinear_local(X0,omega_1,omega_2,omega_3,omega_4)
            % APPLY MAPPING FROM SIM STATES TO MODEL STATES
            % X = [x y z phi theta psi u v w phi_dot theta_dot psi_dot];   % SIM STATE VECTOR
            % Y = [PHI DPHI x Dx y Dy THETA DTHETA PSI DPSI z Dz];         % ODE STATE VECTOR

            % POSITIONS
            % x   = X0(1);     y = X0(2);   z = X0(3);
            PHI = X0(4); THETA = X0(5); % PSI = X(9);
            % VELOCITIES
            Dx   = X0(7); Dy     = X0(8);    Dz = X0(9);
            DPHI = X0(10);DTHETA = X0(11); DPSI = X0(12);
            
            % LOCAL NONLINEAR EXPRESSIONS
            DPSI    = DPSI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            DDPSI   = (27848632988697*omega_2^2)/66701046345491968 - (27848632988697*omega_1^2)/66701046345491968 - (27848632988697*omega_3^2)/66701046345491968 + (27848632988697*omega_4^2)/66701046345491968 + (63122452377224871936*DTHETA*DPHI)/63102186178901703125 + (63122452377224871936*sin(PHI)*tan(THETA)*DTHETA^2)/63102186178901703125 - (63122452377224871936*cos(PHI)*DTHETA*DPHI)/63102186178901703125 + (63122452377224871936*sin(PHI)*DPSI*DPHI)/63102186178901703125 + (63122452377224871936*cos(PHI)*tan(THETA)*DPSI*DTHETA)/63102186178901703125;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
            Dx      = Dx;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
            DDx     = DTHETA*Dz - DPSI*Dy - (300978377846937375*sin(THETA))/30680772461461504;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
            Dy      = Dy;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
            DDy     = (300978377846937375*cos(THETA)*sin(PHI))/30680772461461504 + DPSI*Dx - Dz*DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            DTHETA  = DTHETA;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
            DDTHETA = -(125592194882400484375*omega_2^2*cos(THETA) - 125592194882400484375*omega_1^2*cos(THETA) + 125592194882400484375*omega_3^2*cos(THETA) - 125592194882400484375*omega_4^2*cos(THETA) + 3342365558715433500000*omega_1*cos(THETA)*DPHI - 3342365558715433500000*omega_2*cos(THETA)*DPHI + 3342365558715433500000*omega_3*cos(THETA)*DPHI - 3342365558715433500000*omega_4*cos(THETA)*DPHI + 1033866218355125504000000*cos(THETA)*DPSI*DPHI - 1034198259748452301799424*cos(PHI)*DPSI*DPHI - 1034198259748452301799424*sin(PHI)*DTHETA*DPHI + 1033866218355125504000000*cos(THETA)*cos(PHI)*tan(THETA)*DPSI^2 + 1033866218355125504000000*cos(THETA)*sin(PHI)*tan(THETA)*DPSI*DTHETA + 3342365558715433500000*omega_1*cos(THETA)*cos(PHI)*tan(THETA)*DPSI - 3342365558715433500000*omega_2*cos(THETA)*cos(PHI)*tan(THETA)*DPSI + 3342365558715433500000*omega_3*cos(THETA)*cos(PHI)*tan(THETA)*DPSI - 3342365558715433500000*omega_4*cos(THETA)*cos(PHI)*tan(THETA)*DPSI + 3342365558715433500000*omega_1*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA - 3342365558715433500000*omega_2*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA + 3342365558715433500000*omega_3*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA - 3342365558715433500000*omega_4*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA)/(1034198259748452301799424*cos(THETA));
            Dz      = Dz;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
            DDz     = (300978377846937375*cos(THETA)*cos(PHI))/30680772461461504 - (34446988714643375*omega_1^2)/1005347552017170563072 - (34446988714643375*omega_2^2)/1005347552017170563072 - (34446988714643375*omega_3^2)/1005347552017170563072 - (34446988714643375*omega_4^2)/1005347552017170563072 - Dx*DTHETA + Dy*DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
            DPHI    = DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            DDPHI   = -(125592194882400484375*omega_1^2*cos(THETA) + 125592194882400484375*omega_2^2*cos(THETA) - 125592194882400484375*omega_3^2*cos(THETA) - 125592194882400484375*omega_4^2*cos(THETA) + 1034198259748452301799424*sin(PHI)*DTHETA^2 + 1033866218355125504000000*cos(THETA)*sin(PHI)*DPSI^2 + 1034198259748452301799424*cos(PHI)*DPSI*DTHETA - 1033866218355125504000000*cos(THETA)*cos(PHI)*DPSI*DTHETA - 3342365558715433500000*omega_1*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_2*cos(THETA)*cos(PHI)*DTHETA - 3342365558715433500000*omega_3*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_4*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_1*cos(THETA)*sin(PHI)*DPSI - 3342365558715433500000*omega_2*cos(THETA)*sin(PHI)*DPSI + 3342365558715433500000*omega_3*cos(THETA)*sin(PHI)*DPSI - 3342365558715433500000*omega_4*cos(THETA)*sin(PHI)*DPSI)/(1034198259748452301799424*cos(THETA));                                                                                                                                                                                                                                                                                                                                                                    


%             % LOCAL NONLINEAR EXPRESSIONS
%             DPSI = DPSI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
%             DDPSI = (99326544526713*omega_2^2)/23790539496531427328 - (99326544526713*omega_1^2)/23790539496531427328 - (99326544526713*omega_3^2)/23790539496531427328 + (99326544526713*omega_4^2)/23790539496531427328 + (63122452377224871936*DTHETA*DPHI)/63102186178901703125 + (63122452377224871936*sin(PHI)*tan(THETA)*DTHETA^2)/63102186178901703125 - (63122452377224871936*cos(PHI)*DTHETA*DPHI)/63102186178901703125 + (63122452377224871936*sin(PHI)*DPSI*DPHI)/63102186178901703125 + (63122452377224871936*cos(PHI)*tan(THETA)*DPSI*DTHETA)/63102186178901703125;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
%             Dx = Dx;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%             DDx = DTHETA*Dz - DPSI*Dy - (300978377846937375*sin(THETA))/30680772461461504;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
%             Dy = Dy;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%             DDy = (300978377846937375*cos(THETA)*sin(PHI))/30680772461461504 + DPSI*Dx - Dz*DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
%             DTHETA = DTHETA;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
%             DDTHETA = -(123071212670336359375*omega_2^2*cos(THETA) - 123071212670336359375*omega_1^2*cos(THETA) + 123071212670336359375*omega_3^2*cos(THETA) - 123071212670336359375*omega_4^2*cos(THETA) + 3342365558715433500000*omega_1*cos(THETA)*DPHI - 3342365558715433500000*omega_2*cos(THETA)*DPHI + 3342365558715433500000*omega_3*cos(THETA)*DPHI - 3342365558715433500000*omega_4*cos(THETA)*DPHI + 1033866218355125504000000*cos(THETA)*DPSI*DPHI - 1034198259748452301799424*cos(PHI)*DPSI*DPHI - 1034198259748452301799424*sin(PHI)*DTHETA*DPHI + 1033866218355125504000000*cos(THETA)*cos(PHI)*tan(THETA)*DPSI^2 + 1033866218355125504000000*cos(THETA)*sin(PHI)*tan(THETA)*DPSI*DTHETA + 3342365558715433500000*omega_1*cos(THETA)*cos(PHI)*tan(THETA)*DPSI - 3342365558715433500000*omega_2*cos(THETA)*cos(PHI)*tan(THETA)*DPSI + 3342365558715433500000*omega_3*cos(THETA)*cos(PHI)*tan(THETA)*DPSI - 3342365558715433500000*omega_4*cos(THETA)*cos(PHI)*tan(THETA)*DPSI + 3342365558715433500000*omega_1*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA - 3342365558715433500000*omega_2*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA + 3342365558715433500000*omega_3*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA - 3342365558715433500000*omega_4*cos(THETA)*sin(PHI)*tan(THETA)*DTHETA)/(1034198259748452301799424*cos(THETA));
%             Dz = Dz;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%             DDz = (300978377846937375*cos(THETA)*cos(PHI))/30680772461461504 - (1080177360492106875*omega_1^2)/32171121664549458018304 - (1080177360492106875*omega_2^2)/32171121664549458018304 - (1080177360492106875*omega_3^2)/32171121664549458018304 - (1080177360492106875*omega_4^2)/32171121664549458018304 - Dx*DTHETA + Dy*DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
%             DPHI  = DPHI;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
%             DDPHI = -(123071212670336359375*omega_1^2*cos(THETA) + 123071212670336359375*omega_2^2*cos(THETA) - 123071212670336359375*omega_3^2*cos(THETA) - 123071212670336359375*omega_4^2*cos(THETA) + 1034198259748452301799424*sin(PHI)*DTHETA^2 + 1033866218355125504000000*cos(THETA)*sin(PHI)*DPSI^2 + 1034198259748452301799424*cos(PHI)*DPSI*DTHETA - 1033866218355125504000000*cos(THETA)*cos(PHI)*DPSI*DTHETA - 3342365558715433500000*omega_1*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_2*cos(THETA)*cos(PHI)*DTHETA - 3342365558715433500000*omega_3*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_4*cos(THETA)*cos(PHI)*DTHETA + 3342365558715433500000*omega_1*cos(THETA)*sin(PHI)*DPSI - 3342365558715433500000*omega_2*cos(THETA)*sin(PHI)*DPSI + 3342365558715433500000*omega_3*cos(THETA)*sin(PHI)*DPSI - 3342365558715433500000*omega_4*cos(THETA)*sin(PHI)*DPSI)/(1034198259748452301799424*cos(THETA));                                                                                                                                                                                                                                                                                                                                                                    

            % APPLY MAPPING FROM THE MODEL STATES TO THE SIM STATES
            dX = [Dx Dy Dz DPHI DTHETA DPSI DDx DDy DDz DDPHI DDTHETA DDPSI]';
        end
        % GET BODY PARAMETERS
        function [BODY] = getBodyParameters(config)
            % This function gets the phyiscal parameters of the ARdrone
            % airframe, as measured from the physical system.
            
            switch upper(config)
                case 'NONE'
                    BODY = struct('mass',0.366,'width',0.450,'length',0.290);
                case 'OUTDOOR'
                    BODY = struct('mass',0.400,'width',0.452,'length',0.452);
                case 'INDOOR'
                    BODY = struct('mass',0.436,'width',0.515,'length',0.515);
                otherwise
                    error('Configuration not recognised');
            end
            % INERTIAL MATRIX (x) configuration (FICTIONAL I MATRIX)
            BODY.I = 1.0E-06.*[2.8032E+04,0.0000E+00,0.0000E+00;
                               0.0000E+00,2.8032E+04,0.0000E+00;
                               0.0000E+00,0.0000E+00,2.8023E+04];  
        end
    end
end