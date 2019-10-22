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
        function [obj] = agent_2D(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                            % Call the super class            

            % Assign defaults (2D)
            obj = obj.SetDetectionRadius(inf);
            obj = obj.SetVIRTUALparameter('is3D',false);    % Indicate the agent is operating in 2d
            obj.DYNAMICS = obj.CreateDYNAMICS();            % Override 3D dynamics with 2D
            obj.localState = zeros(6,1);                    % Default state
            
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin);             % Recursive overrides
            % Re-affirm associated properties
            obj.MEMORY = obj.GetMemoryStructure(obj.maxSamples);  % Initialise memory structure (with 2D varient)             
            % /////////////////////////////////////////////////////////////
        end
        % Setup (state [x y psi])
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.

            % Initialises the state vector [x y psi; dx dy dpsi]
            [obj] = obj.initialise_2DVelocities(localXYZVelocity,localXYZrotations);
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
                                    
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(ENV,varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredHeadingVector = obj.GetTargetHeading();                 % Design the current desired trajectory from the waypoint.  
            desiredVelocity = desiredHeadingVector*obj.nominalSpeed;
                  
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
%             [obj] = obj.updateGlobalProperties(ENV.dt,newState);

            [obj] = obj.controller(ENV.dt,desiredVelocity);
        end
    end
    
    %% /////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods   
        % CONTROLLER (2D)
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired 2D velocity vector. This function is designed
            % to act as a 2D override for the 3D controller function be the
            % same name in 'agent'.
            
            % Input sanity check #1 - Is feasible
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 2,'Requested velocity vector is not a 2D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.'); 
            
            % Input sanity check #2 - Zero vector
            [heading,speed] = obj.nullVelocityCheck(desiredVelocity);   

            % Get the relative heading
            [dHeading] = obj.GetVectorHeadingAngles([1;0],heading); % Relative heading angles   
            % The desired heading rate
            omega = (-dHeading)/dt;
            
            % Apply kinematic constraints to state changess
            [omega_actual,speed_actual] = obj.ApplyKinematicContraints(dt,omega,speed);
            
            % OMIT TRAJECTORY CHANGES IF IDLE
            if obj.GetIdleStatus()
                omega_actual = 0;
                speed_actual = 0;
                obj = obj.SetNominalSpeed(0);
            end
            
            % //////////////////// SIMPLE DYNAMICS ////////////////////////
            [dX] = obj.dynamics_simple(obj.localState(1:3),[speed_actual;0],omega_actual);
            obj.localState(1:3) = obj.localState(1:3) + dt*dX;
            obj.localState(4:6) = dX; 
                        
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES ///////////////
            [obj] = obj.updateGlobalProperties_2DVelocities(dt,obj.localState);
        end
        % IMPOSE 2D KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = kinematicContraints(obj,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % HEADING IS RELATIVE
            currentVelocity = obj.localState(4:5,1);
            currentHeading = zeros(size(desiredHeading));
            
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - currentHeading;
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = obj.boundValue(Hdot,-obj.DYNAMICS.angularVelocityLimits,obj.DYNAMICS.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = currentHeading + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed;                          % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = obj.boundValue(acceleration,-obj.DYNAMICS.linearAccelerationLimits(1),obj.DYNAMICS.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = obj.boundValue(achievedSpeed,-obj.DYNAMICS.linearVelocityLimits(1),obj.DYNAMICS.linearVelocityLimits(1)); 
        end
    end
    %% ////////////////// GLOBAL UPDATE/STATE FUNCTIONS ///////////////////
    methods
        % INITIAL 2D STATE - [x y psi dx dy dpsi]
        function [obj] = initialise_2DVelocities(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            
            % BUILD STATE VECTOR
            obj.localState = zeros(6,1);
            obj.localState(3) = localXYZrotations(3);
            obj.localState(4:5) = localXYZVelocity(1:2);
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj = obj.SetVIRTUALparameter('priorState',obj.localState);
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x y psi dx dy dpsi]'
        function [obj] = updateGlobalProperties_2DVelocities(obj,dt,eulerState)
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
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = obj.GetVIRTUALparameter('quaternion');   
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(quaternion_k');
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*localAngularRates;
            % UPDATE THE QUATERNION POSE
            quaternion_k_plus = OMAS_geometry.integrateQuaternion(quaternion_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(quaternion_k_plus');
                        
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*localLinearRates;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end
    end
end
