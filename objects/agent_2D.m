%% THE 2D AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the 'agent' class, this object class is intended provide
% support for 2D simulation. The functions provided here 

% Author: James A. Douthwaite

classdef agent_2D < agent
%%% AGENT(2D) BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % PROPERTIES UNIQUE TO THE RECIPROCAL VO METHOD
    end
%%   CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function [obj] = agent_2D(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                        % Call the super class
            % AGENT SENSES IN 3D
            obj.VIRTUAL.is3D = logical(false);          % Indicate the agent is operating in 2d
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % ///////////////////// SETUP FUNCTION ////////////////////////////
        % DEFAULT 2D STATE - [x y psi]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi]
            
            % BUILD STATE VECTOR
            obj.localState = zeros(3,1);
            obj.localState(3) = localXYZrotations(3);
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj.VIRTUAL.priorState = obj.localState;
        end
        % ///////////////// AGENT MAIN CYCLE //////////////////////////////
        function [obj] = main(obj,TIME,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % DEFAULT BEHAVIOUR 
            desiredSpeed = 1;
            desiredHeadingVector = [1;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
                desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
            end
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            headingRate = 0.0;      % Constant +ve yaw rate
            speed = 1;              % Constant speed (m/s)
            
            if TIME.currentTime >= 0 && TIME.currentTime < 5
                headingRate = 0.5;
                fprintf('right...\n');
            end
            
            if TIME.currentTime >= 5 && TIME.currentTime < 10
                headingRate = -0.5;
                fprintf('left...\n');
            end
            
            if TIME.currentTime >= 10 && TIME.currentTime < 15
                headingRate = 0.5;
                fprintf('right...\n');
            end   
            
            % GET THE NEW STATE VECTOR
            [dXdt] = obj.dynamics_singleIntegrator(obj.localState,[speed;0],headingRate);
            newState = obj.localState + dXdt*dt;
            % UPDATE THE CLASS GLOBAL PROPERTIES 
            [obj] = obj.updateGlobalProperties_ENU(dt,newState);
        end
    end
    
    % ///////////////////// (2D) SENSORY FUNCTIONS ////////////////////////
    methods
        % 2D SENSOR MODEL - LOCAL GPS & PITOT TUBE
        function [position,velocity,radius] = getAgentMeasurements(obj)
            % This function makes estimates of the agent's current position
            % and velocity states (defined in obj.localState).
            positionIndices = 1:2;
            velocityIndices = 4:5;
            % GET THE LOCAL ABSOLUTE VARIABLES
            position =  obj.localState(positionIndices,1) + obj.SENSORS.sigma_position*rand(2,1); % Get the absolute position measurement
            velocity =  obj.localState(velocityIndices,1) + obj.SENSORS.sigma_velocity*rand(2,1); % Get the absolute velocity measurement          
            % RADIUS IS ASSUMED KNOWN
            radius = obj.VIRTUAL.radius;
        end
    end
    
    % //////////////////////// DYNAMICS & CONTROL //////////////////////////
    methods   
        % 2D CONTROLLER
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % NUMERICAL SANITY CHECK ONE
            assert(~any(isnan(desiredVelocity)) && isnumeric(desiredVelocity),'Requested velocity vector must be a 2D local vector');
            % NUMERICAL SANITY CHECK TWO
            [unitDirection,desiredSpeed] = obj.nullVelocityCheck(desiredVelocity);   
            
            % APPLY SPEED CONSTRAINT
            if abs(desiredSpeed) > obj.maxSpeed
                desiredSpeed = sign(desiredSpeed)*obj.maxSpeed; 
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dHeading] = obj.getVectorHeadingAngles([1;0],unitDirection); % Relative heading angles   
            omega = -dHeading/dt;
            desiredVelocity = [desiredSpeed;0];
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:3),desiredVelocity,omega);
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
        % INITIALISE 2D DYNAMIC PARAMETERS
        function [obj] = getDynamicParameters(obj)
            % This function calls the dynamics structure that gives the
            % agent dynamic/kinematic limits and parameters.
            % BUILD THE DYNAMICS SUB-STRUCTURE
            obj.DYNAMICS = struct('linearVelocityLimits',[inf;inf],...     % Limits on the agents linear velocity 
                              'linearAccelerationLimits',[inf;inf],...     % Limits on the agents linear acceleration
                                 'angularVelocityLimits',[inf;inf],...     % Limits on the agents angular velocity
                             'angularAccelerationLimits',[inf;inf]);       % Limits on the agents angular acceleration
        end
    end
    
    % /////////// AGENT/DERIVATES GLOBAL UPDATE/STATE FUNCTIONS ///////////
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
            obj.VIRTUAL.priorState = obj.localState;
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x y psi dx dy dpsi]'
        function [obj] = updateGlobalProperties_2DVelocities(obj,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y psi dx dy dpsi].
 
            velocityIndices = 4:5;
            headingRateIndex = 6;
            
            % EQUIVALENT RATES
            velocity_k_plus = eulerState(velocityIndices);                 % Use the linear velocity directly.
            localAxisRates  = eulerState(headingRateIndex);                % Use the heading rate directly.
                        
            % MAP TO GLOBAL 3D REPRESENTATION
            velocity_k_plus = [velocity_k_plus;0];  
            localAxisRates  = [0;0;1]*localAxisRates;
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = obj.VIRTUAL.quaternion;   
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(quaternion_k');
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*localAxisRates;
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
                                                      eulerState);              % The new state for reference
        end
        % [DEFINED IN 'agent.m']
        % GLOBAL UPDATE - EULER 6DOF(3DOF) [x y z phi theta psi],[x y psi]
        %function: updateGlobalProperties_ENU(obj,dt,eulerState)
        % GLOBAL UPDATE - EULER 6DOF(3DOF) (NO ROTATIONS)
        %function: updateGlobalProperties_ENU_fixedFrame(obj,dt,eulerState)
    end
end
