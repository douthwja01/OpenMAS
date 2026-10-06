%% THE 2D AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the 'agent' class, this object class is intended provide
% support for 2D simulation. The functions provided here 

% Author: James A. Douthwaite

classdef agent_2D < agent
    % AGENT(2D) BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.

    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_2D(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            this@agent(varargin);                    % Call the super class            

            % Assign required 2D overrides
            this.localState = zeros(6,1);                       % Default state
            this.SetGLOBAL('is3D',false);                       % Indicate the agent is operating in 2d
            this.DYNAMICS = this.CreateDYNAMICS();              % Override 3D dynamics with 2D
            this.MEMORY	  = this.CreateMEMORY(this.maxSamples); % Override 3D memory structure with a 2D memory varient

            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides         
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        % - The same as any 2D/3D object
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
                                    
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredHeadingVector = this.GetTargetHeading();                 % Design the current desired trajectory from the waypoint.  
            desiredVelocity = desiredHeadingVector*this.v_nominal;
                  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            headingRate = 0.0;      % Constant +ve yaw rate
            speed = 1;              % Constant speed (m/s)
            
            if ENV.currentTime >= 0 && ENV.currentTime < 5
                headingRate = 0.5;
                fprintf('right...\n');
            end
            
            if ENV.currentTime >= 5 && ENV.currentTime < 10
                headingRate = 0.5;
                fprintf('left...\n');
            end
            
            if ENV.currentTime >= 10 && ENV.currentTime < 15
                headingRate = 0.5;
                fprintf('right...\n');
            end   
            
            if ENV.currentTime == 1
                
            end
            
%             % GET THE NEW STATE VECTOR
%             [dXdt] = obj.dynamics_singleIntegrator(obj.localState,[speed;0],headingRate);
%             newState = obj.localState + dXdt*ENV.dt;
%             % UPDATE THE CLASS GLOBAL PROPERTIES 
%             [obj] = obj.UpdateGlobalProperties(ENV.dt,newState);

            [this] = this.Controller(ENV.dt,desiredVelocity);
        end
    end
    
    %% /////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods   
        % CONTROLLER (2D)
        function [this] = Controller(this,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired 2D velocity vector. This function is designed
            % to act as a 2D override for the 3D controller function be the
            % same name in 'agent'.
            
            % Input sanity check #1 - Is feasible
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 2,'Requested velocity vector is not a 2D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.'); 
            
            % Input sanity check #2 - Zero vector
            [heading,speed] = this.nullVelocityCheck(desiredVelocity);   

            % Get the relative heading
            [dHeading] = this.GetVectorHeadingAngles([1;0],heading); % Relative heading angles   
            % The desired heading rate
            omega = (-dHeading)/dt;
            
            % Apply kinematic constraints to state changess
            [omega_actual,speed_actual] = this.ApplyKinematicContraints(dt,omega,speed);
%             omega_actual = omega;
%             speed_actual = speed;
            
            % OMIT TRAJECTORY CHANGES IF IDLE
            if this.IsIdle()
                omega_actual = 0;
                speed_actual = 0;
                this.v_nominal = 0;
            end
            
            % //////////////////// SIMPLE DYNAMICS ////////////////////////
            [dX] = this.SingleIntegratorDynamics(this.localState(1:3),[speed_actual;0;omega_actual]);
            this.localState(1:3) = this.localState(1:3) + dt*dX;
            this.localState(4:6) = dX; 
                        
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES ///////////////
            [this] = this.GlobalUpdate_2DVelocities(dt,this.localState);
        end
        % IMPOSE 2D KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = KinematicContraints(this,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % HEADING IS RELATIVE
            currentVelocity = this.localState(4:5,1);
            currentHeading = zeros(size(desiredHeading));
            
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - currentHeading;
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = this.boundValue(Hdot,-this.DYNAMICS.angularVelocityLimits,this.DYNAMICS.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = currentHeading + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed;                          % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = this.boundValue(acceleration,-this.DYNAMICS.linearAccelerationLimits(1),this.DYNAMICS.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = this.boundValue(achievedSpeed,-this.DYNAMICS.linearVelocityLimits(1),this.DYNAMICS.linearVelocityLimits(1)); 
        end
    end
    %% ////////////////// GLOBAL UPDATE/STATE FUNCTIONS ///////////////////
    methods
        % INITIAL 2D STATE - [x y psi dx dy dpsi]
        function [this] = setup_2DVelocities(this,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            
            % BUILD STATE VECTOR
            this.localState      = zeros(6,1);
            this.localState(3)   = localXYZrotations(3);
            this.localState(4:5) = localXYZVelocity(1:2);
            % RETAIN THE PRIOR STATE FOR REFERENCE
            this.SetGLOBAL('priorState',this.localState);
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x y psi dx dy dpsi]'
        function [this] = GlobalUpdate_2DVelocities(this,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y psi dx dy dpsi].
 
            velocityIndices  = 4:5;
            headingRateIndex = 6;
            
            % EQUIVALENT RATES
            localLinearRates  = eulerState(velocityIndices);               % Use the linear velocity directly.
            localAngularRates = eulerState(headingRateIndex);              % Use the heading rate directly.
                        
            % MAP TO GLOBAL 3D REPRESENTATION
            localLinearRates  = [localLinearRates;0];  
            localAngularRates = [0;0;1]*localAngularRates;
            % Get the current parameter
            p_k = this.GetGLOBAL('position');  
            v_k = this.GetGLOBAL('velocity');  
            q_k = this.GetGLOBAL('quaternion');  
            
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
             
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(q_k');
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*localAngularRates;
            % UPDATE THE QUATERNION POSE
            q_k_plus = OMAS_geometry.integrateQuaternion(q_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(q_k_plus');
                        
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            v_k_plus = R_k_plus'*localLinearRates;
            p_k_plus = p_k + dt*v_k;
            
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.GlobalUpdate_direct(...
                p_k_plus,...    % Global position at k plius
                v_k_plus,...    % Global velocity at k plus
                q_k_plus);      % Quaternion at k plus
        end
    end
end
