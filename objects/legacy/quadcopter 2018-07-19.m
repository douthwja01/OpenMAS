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
    %  CLASS METHODS
    methods
        % CONSTRUCTOR METHOD
        function obj = quadcopter(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                         % Create the super class 'agent'
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            % IMPORT THE AGENT DYNAMICS & CONTROL PARAMETERS
            [obj.DYNAMICS] = obj.importModelProperties();
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
        function [obj] = processTimeCycle(obj,TIME,varargin)
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
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,~] = obj.getAgentUpdate(varargin{1});                   % IDEAL INFORMATION UPDATE 
         
            % DETERMINE THE STATE REFERENCE BASED ON CONTROL METHOD
            [obj,X_ref] = obj.setControlType('position',obj.localState);
            
            % //////// ///// UPDATE THE LOCAL STATE VECTOR ////////////////
            % Continous feedback to a desired state
            [obj.localState] = obj.updateLocalState(TIME,obj.localState,X_ref);
            
            % GET THE LATEST ROTATION MATRIX
            Rnew = eul2rotm(obj.localState(7:9)','XYZ');                      
            quaternion = OMAS_axisTools.rotationMatrixToQuaternion(Rnew');

            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            [obj] = obj.updateGlobalProperties_direct(obj.localState(1:3),...
                                                      obj.localState(4:6),...
                                                      quaternion,...
                                                      obj.localState);   
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputNames = {'p (rad/s)','q (rad/s)','r (rad/s)'};
            obj.DATA.inputs(1:3,TIME.currentStep) = obj.localState((end-2):end);   % Record the control inputs
        end
    end
    methods
        % CONTROL METHODOLOGY
        function [obj,X_ref] = setControlType(obj,type,X)
            % DEFAULT BEHAVIOUR
            p = zeros(3,1);
            v = zeros(3,1);
            eta = zeros(3,1);
            omega = zeros(3,1);
            
            % CURRENT ROTATIONS
            R = OMAS_axisTools.eulersToRotationMatrix(X(7:9));
            
            % WAYPOINT BEHAVIOUR
            if ~isempty(obj.targetWaypoint)
                desiredSpeed = 2;
            	p_wp = obj.targetWaypoint.position;                        % Target relative position 
                desiredHeadingVector = OMAS_axisTools.unit(p_wp);
            else
                desiredSpeed = 0;
                p_wp = zeros(3,1);                                         % Default target relative position
                desiredHeadingVector = [1;0;0];
            end
            
            % DEFINE THE CONTROLLER METHODOLOGY
            switch type
                case 'position'
                    % DIRECTLY ALLOCATE THE TARGET WAYPOINT POSITION
                    p = R*p_wp + X(1:3); % Relative position + global position
                case 'velocity' 
                    % REMOVE POSITIONAL FEEDBACK
                    obj.DYNAMICS.K(1:4,1:3) = zeros(4,3); 
                    % DEFINE THE DESIRED GLOBAL VELOCITY
                    v = R'*(desiredHeadingVector*desiredSpeed);
                case 'heading'
                    % REMOVE POSITIONAL FEEDBACK
                    obj.DYNAMICS.K(1:4,1:3) = zeros(4,3); 
                    % DEFINE THE GLOBAL VELOCITY 
                    v = R'*([1;0;0]*desiredSpeed);
                    % GLOBAL HEADING ANGLE
                    
                    
                    eta = [0;0;0];
                otherwise
                    error('Control methodology not recognised');
            end
            % CREATE THE REFERENCE STATE 
            X_ref = [p;v;eta;omega]; 
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
                [~,Xset] = ode45(@(t,X) obj.stateDynamics(X,X_desired), tspan, X0, opts);
                X = Xset(end,:)';
            end
        end    
        % CONTINUOUS FEEDBACK - CLOSED LOOP
        function [dX] = stateDynamics(obj,X,X_desired)
            % PARSE THE STATE VECTOR
            % X = [x y z phi theta psi dx dy dz dphi dtheta dpsi]
            pn = X(1:3);
            vn = X(4:6);
            en = X(7:9);
            omega = X(10:12);
            
            % DYNAMICS PROPERTIES
            g  = obj.DYNAMICS.g;
            e3 = obj.DYNAMICS.e3;
            m  = obj.DYNAMICS.m;
            M  = obj.DYNAMICS.M;
            K  = obj.DYNAMICS.K;
               
            R = OMAS_axisTools.eulersToRotationMatrix(en);
            
            % STATE ERROR
            delta_x = [pn;vn;en;omega] - X_desired;

            % INPUTS
            inputs = [m*g;0;0;0] - K*delta_x;
            
            % INPUTS ARE FORCES AND TORQUES DIRECTLY
            f=inputs(1);
            tau=inputs(2:4);
            
           	% NON-LINEAR DYNAMICS
            % TRANSLATIONAL
            pn_dot = vn;
            vn_dot = f/m*R*e3-g*e3;
            % ROTATIONAL
            en_dot = R*omega;
            omega_dot = inv(M)*(tau - OMAS_axisTools.skew(omega)*M*omega);
            % REBUILD THE STATE DIFFERENTIAL
            dX=[pn_dot;vn_dot;en_dot;omega_dot];
        end
        % INTIALISE - QUADCOPTER STATE VECTOR [pn;vn;etn;omega]
        function [obj] = initialise_localState(obj,localXYZVelocity,localXYZrotations)
            % This function is called in order to build a local state
            % vector intended for aerospace control simulation (NED).
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z dx dy dz R1-R9 p q r]

            % THIS MODEL OPERATES IN THE GLOBAL FRAME
            % Properties can be taken directly from the VIRTUAL properties
            pn0 = obj.VIRTUAL.globalPosition;       % True global position
            % vn0 = obj.VIRTUAL.globalVelocity;       % True global velocity
            % qn0 = obj.VIRTUAL.quaternion;           % True global quaterion          
            
            % GET THE ROTATION MATRIX fixed local axes -> rotated global pose
            % R0eval = OMAS_axisTools.quaternionToRotationMatrix(qn0); % For validation
            R0 = OMAS_axisTools.eulersToRotationMatrix(localXYZrotations);
            
            % BUILD THE GLOBAL STATE VECTOR     
            vn0 = R0'*localXYZVelocity;               % The global velocity
            omega0 = zeros(3,1);
            x0 = [pn0;vn0;-localXYZrotations;omega0];                      
            % ASSIGN THE LOCAL FRD STATE
            obj.VIRTUAL.priorState = x0;
            obj.localState = x0;
        end
    end
    % STATIC FUNCTIONS
    methods (Static)
        % IMPORT THE MODEL PROPERTIES
        function [modelParams] = importModelProperties()
            
            % BUILD THE DYNAMIC PROPERTIES OF THE AGENT
            g=9.8;
            I=eye(3);
            e3=I(:,3);
            m=2*0.5;%kg
            M=2*diag([4.856, 4.856, 8.801])*10^(-3);
            psi_desired = 0;
            
            % PLANT MATRIX
            A=zeros(12);
            A(1:3,4:6) = eye(3);
            A(4:6,7:9) = g*[sin(psi_desired), cos(psi_desired), 0;
                           -cos(psi_desired), sin(psi_desired), 0;
                                           0,                0, 0];
            A(7:9,10:12) = eye(3);
            % INPUT MATRIX
            B=zeros(12,4); % input matrix
            B(4:6,1)=e3/m;
            B(10:12,2:4)=inv(M);
            % OTHER STATE SPACE MATRICES
            C = eye(12);
            D = [zeros(8,4);eye(4)];
            % LQR CONTROLLER                       
            isdescrete = 0;
            if isdescrete
                SYS = ss(A,B,C,D);
                [SYSd] = c2d(SYS,dt,'zoh');
                Q = eye(12)*1E0;
                R = eye(4)*1E4;  
                % LQR DESCRETE FEEDBACK
                [K_lqr] = dlqr(SYSd.A,SYSd.B,Q,R);
            else
                Q = eye(12);
                R = eye(4);
                % LQR CONTINUOUS FEEDBACK
                [K_lqr,~,~] = lqr(A,B,Q,R);
            end
            
            % RE-ORGANISE THE DYNAMIC PARAMETERS
            modelParams = struct('g',g,...      % Gravitational constant
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