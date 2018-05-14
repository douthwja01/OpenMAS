%% GENERIC OBJECT CLASS (objectDefinition.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic object and import its variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite 20/05/2016

classdef objectDefinition
%%% BASIC OBJECT CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % SIMULATION IDENTIFIERS
        objectID;
        name;
        % KINEMATIC PROPERTIES
        localState  = zeros(12,1);   % Other properties derived from the state.
        % VIRTUAL PROPERTIES
        VIRTUAL;                     % Virtual (SIMULATION) data container
    end
%   CLASS METHODS
    methods (Access = public)
        % CONSTRUCTOR METHOD
        function obj = objectDefinition(varargin)
            % This function is designed to automatically apply a generic
            % naming convention to all objects generated prior to a
            % simulation.
            % INPUTS:
            % varargin - A cell array of property/value pairs 
            % OUTPUTS:
            % obj     - The generated object
            
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            
            % GENERATE OBJECT-ID (If NOT ALLOCATED BY A SUPERCLASS)
            spawnNewObject = 0;
            if isempty(obj.objectID)
                spawnNewObject = 1;
                [obj.objectID] = objectDefinition.getObjectID();           % Assign the number as an ID tag
            end
            
            % GENERATE OBJECT NAME
            if isempty(obj.name) || strcmp(obj.name,'')
                % No input name string specified, use greek naming scheme
                obj.name = obj.getObjectName(obj.objectID);
            end
                        
            % DEFINE DEFAULT VIRTUAL PROPERTIES
            obj.VIRTUAL.type = OMAS_objectType.misc;                       % Indicate whether the agent is 'active' or 'passive'
            obj.VIRTUAL.radius = 1;                                        % Virtual size
            obj.VIRTUAL.colour = rand(1,3,'single');                       % Virtual colour (for plotting)
            obj.VIRTUAL.symbol = 'square';
            obj.VIRTUAL.globalPosition = [0;0;0];
            obj.VIRTUAL.globalVelocity = [0;0;0];
            obj.VIRTUAL.quaternion = [1;0;0;0];
            obj.VIRTUAL.idleStatus = logical(true);
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
            
            % DISPLAY THE CURRENT AGENT PROPERTIES
            if spawnNewObject
                fprintf('objectcount: %d\tname: %s\ttype: %s\n',obj.objectID,obj.name,class(obj));
            end
        end 
        % ////////////////// THE DEFAULT OBJECT CYCLE /////////////////////
        function obj = processTimeCycle(obj,TIME,varargin)
            % This is a generic process cycle of an object that accepts no
            % input commands/feedback and simply updates its states based
            % on its current attributes
            % INPUT:
            % varargin - Generic variable input container, containing the simulation timestep
            % OUTPUT:
            % obj - The updated object
            
            % DETERMINE THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % MAINTAIN LOCAL INPUT CONDITIONS
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_singleIntegrator(dt,[0;0;0],[0;0;0]);
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
    end
    % /////////////////////// DYNAMICS & CONTROL /////////////////////////
    methods
        % The agent object has both a local NED and global ENU state to
        % allow control and avoidance simulation to be conducted
        % simulaneously within the common frame of references. The local
        % state vector is of the form:
        % Xi_frd = [x;y;z;psi;the;phi;u;v;w;p;q;r]_i
        % where as the global is of the form:
        % Xi_enu = [x;y;z;dx;dy;dz;q1;q2;q3;q4]_i.
        
        % The functions below provide the utilities to update both during
        % each time cycle.
        
        % UPDATE AGENT STATE VECTOR AS A TWO WHEELED ROBOT
        function state_k_plus = stateDynamics_differencialDrive(obj,dt,a_left,a_right)
            % Becker, M. (2006). Obstacle avoidance procedure for mobile robots. ABCM Symposium Series
            % [x;y;theta;v;omega]
            state_k = obj.localState; % Get the current state
            % Build the state difference
            dX = zeros(5,1);
            dX(1) = state_k(4)*cos(state(3)); % x velocity
            dX(2) = state_k(4)*sin(state(3)); % y velocity
            dX(3) = state_k(5);               % rotational rate
            dX(4) = (a_right + a_left)/2;                   % Linear acceleration
            dX(5) = (a_right - a_left)/obj.wheelSeperation; % Angular acceleration
            % GENERATE THE NEXT STATE
            state_k_plus = state_k + dt*dX;
        end
        % UPDATE AGENT STATE VECTOR AS A UNI-CYCLE ROBOT
        function state_k_plus = stateDynamics_unicycle(obj,dt,speed,headingRate)
            % This function provides the dynamics of a uni-cycle model
            % moving in on a 2D plane.
            
            % THE PREVIOUS STATE
            state_k = obj.localState;
            % THE STATE UPDATE
            dX = zeros(3,1);
            dX(1) = speed*cos(state_k(3));
            dX(2) = speed*sin(state_k(3));
            dX(3) = headingRate;
            % CALCULATE NEW STATE
            state_k_plus = state_k + dt*[dX;zeros(3,1)]; % Increment the state
            state_k_plus(4:6,1) = dX;                    % retain the velocities
        end
        
        % UPDATE AGENT STATE VECTOR FROM LOCAL ACCELERATIONS
        function state_k_plus = stateDynamics_doubleIntegrator(obj,dt,linearAcceleration,headingAcceleration)
            % This function provides a basic model for the agent
            % kinematics, controlled by linear and angular acceleration inputs
            % INPUTS:
            % dt                  - The simualation timestep
            % linearAcceleration  - The agent's current linear acceleration [2x1] or [3x1]
            % angularAcceleration - The agent's current angular acceleration [1x1] or [3x1]
            % OUTPUTS:
            % newState        - The updated local (FRD) state vector [12x1]
            
            % LOCAL STATE:
            % [x; y; z; phi; theta; psi; u; v; w; phi_dot; theta_dot; psi_dot]
            
            % SET THE INPUT ORDER
            % The input "linearAcceleration" and "headingAcceleration" are
            % related to the second order of the states directly.
            systemOrder = 2;
            
            % CALL THE N-INTEGRATOR
            state_k_plus = obj.stateDynamics_NIntegrator(systemOrder,dt,linearAcceleration,headingAcceleration);
        end
        % UPDATE AGENT STATE VECTOR FROM LOCAL VELOCITIES
        function state_k_plus = stateDynamics_singleIntegrator(obj,dt,linearRates,headingRates)
            % This function preforms the state update directly from a
            % velocity vector in the body axis system.
            % ASSUMPTIONS:
            % - The agent is capable of changing its heading instanaeously.
            % - The heading rates is the euler heading rates.
            
            % INPUTS:
            % linearRates  - Defined as the rates along the body axes (m/s)
            % headingRates - Defined as the euler rates [dphi;dtheta;dpsi] (rad/s)
            
            % SET THE INPUT ORDER
            % The input "linearRates" and "headingRate" are related are the
            % first order of the states directly.
            systemOrder = 1;
            
            % CALL THE N-INTEGRATOR
            state_k_plus = obj.stateDynamics_NIntegrator(systemOrder,dt,linearRates,headingRates);
        end
        % GET UPDATE FOR AN 'N'INTEGRATOR SYSTEM, INPUT OF THAT ORDER
        function state_k_plus = stateDynamics_NIntegrator(obj,systemOrder,dt,nOrderLinearInputs,nOrderAngularInputs)
            % This function preforms the state update from an input vector
            % of relation N to the given state parameter.
            % ASSUMPTIONS:
            % - The agent is capable of changing its heading instanaeously.
            % - The heading rates is the euler heading rates.
            % - obj.localState is of length N*states of the form [x,dx,ddx.]
            
            % INPUTS:
            % systemOrder  - If the system (and input) is zeroth order, first order etc...
            % linearRates  - Defined as the rates along the body axes (m/s)
            % headingRates - Defined as the euler rates [dphi;dtheta;dpsi] (rad/s)
            
            % CHECK EXPECTED INPUT DIMENSIONS
            Aflag = numel(nOrderLinearInputs) == 3 && numel(nOrderAngularInputs) == 3; % IS 3D
            Bflag = numel(nOrderLinearInputs) == 2 && numel(nOrderAngularInputs) == 1; % IS 2D
            assert((Aflag + Bflag) == 1,'Integrator inputs are of the wrong dimensions');
            
            % CALCULATE THE STATE-INPUT PARAMETERS
            % Either rotation occurs in one dimension (2D) or in three (3D)
            rotationDOF = numel(nOrderAngularInputs); % Is heading a yaw (2D) or [phi,theta psi] (3D)
            if rotationDOF == 1
                linearIndices = 1:2;
                angularIndices = 3;
            else
                linearIndices = 1:3;
                angularIndices = 4:6;
            end
            stateNumber = numel([linearIndices,angularIndices]);            % The length of the state vector
            
            % COMPARE TO THE OBJECTS CURRENT STATE
            assert(numel(obj.localState) == stateNumber*(systemOrder+1),...
                'The objects local state vector is of different order to the update.');
            
            % THE STATE DEVIATION (DESCRIBING SYSTEM INPUTS)
            inputs = zeros(stateNumber,1);
            inputs(linearIndices,1)  = nOrderLinearInputs;
            inputs(angularIndices,1) = nOrderAngularInputs;      % The 'nth' order axial inputs
            
            % THE PREVIOUS STATE
            state_k = obj.localState;
            % THE STATE DEVIATION
            dX = [state_k((stateNumber+1):stateNumber*(systemOrder+1));zeros(stateNumber,1)];
            % BUILD THE NEW STATE UPDATE (% X_(k+1) = X_(k) + dX_k)
            state_k_plus = zeros(stateNumber*(systemOrder+1),1);
            state_k_plus(1:stateNumber*(systemOrder+1),1) = state_k + dt*dX;
            state_k_plus(stateNumber*systemOrder+1:end,1) = inputs;
        end
        
        % UPDATE AGENT STATE VECTOR FROM A VELOCITY VECTOR
        function state_k_plus = stateDynamics_simple(obj,dt,velocity_k_plus,omega_k_plus)
            % This function assumes that the velocity changes are
            % implemented this timestep directly, integration then occurs
            % including the updates from this timestep:
            % State vector is assumed:
            % 2D - [x y psi dx dy dpsi]'
            % 3D - [x y z phi theta psi dx dy dz dphi dtheta dpsi]'
            
            % CHECK EXPECTED INPUT DIMENSIONS
            Aflag = numel(velocity_k_plus) == 3 && numel(omega_k_plus) == 3; % IS 3D
            Bflag = numel(velocity_k_plus) == 2 && numel(omega_k_plus) == 1; % IS 2D
            assert((Aflag + Bflag) == 1,'Integrator inputs are of the wrong dimensions');
            
            % Either rotation occurs in one dimension (2D) or in three (3D)
            rotationDOF = numel(omega_k_plus); % Is heading a yaw (2D) or [phi,theta psi] (3D)
            if rotationDOF == 1
                linearIndices = 1:2;
                angularIndices = 3;
            else
                linearIndices = 1:3;
                angularIndices = 4:6;
            end
            stateNumber = numel([linearIndices,angularIndices]);
            
            % CALCULATE THE STATE UPDATE
            state_k_plus = zeros(2*stateNumber,1);
            state_k_plus(linearIndices)  = obj.localState(linearIndices)  + dt*velocity_k_plus;
            state_k_plus(angularIndices) = obj.localState(angularIndices) + dt*omega_k_plus;
            state_k_plus(stateNumber+linearIndices)  = velocity_k_plus;
            state_k_plus(stateNumber+angularIndices) = omega_k_plus; % Record the input velocities
        end
    end
    
    % ////////////////// SIMULATION & CORE INTERFACES /////////////////////
    methods
        % INTIALISE LOCAL 12DOF ENU STATE VECTOR 
        function [obj] = initialise_localState(obj,localENUVelocity,localENUrotations)
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
            
            % BUILD THE DEFAULT ENU STATE VECTOR
            localFLUState = zeros(12,1);
            localFLUState(4:6) = localENUrotations;                        % Get the initial local euler heading
            localFLUState(7:9) = localENUVelocity;                         % Get the initial local velocity vector
            
            % CONVERT THE FLU (LOCAL ENU) STATES TO THE FRD (LOCAL NED)
            %             [obj.localState] = OMAS_axisTools.convertNEDstateToENU(localFLUState); % Assign the FRD(NED) state
            obj.localState = localFLUState;
        end
        % THE GLOBAL PROPERTY FOR GENERAL OBJECT CLASS DERIVATIVES
        function [obj] = updateGlobalProperties(obj,dt,localENUState)
            % This function is used to calculate the new global (.VIRTUAL)
            % parameters for the current object.

            % updates the agents global properties from
            % X_global(@t=k] to X_global(@t=k+1)
            
            % globalPosition - The global ENU position
            % globalVelocity - The global ENU velocity
            % quaternion - Rotation between the local ENU and global ENU
            
            % CONSTANT PARAMETERS
            velocityIndices = 7:9;
            omegaIndices = 10:12;
            
            % DEFINE UPDATE PARAMETERS
%             ENUrates = localENUState(10:12); % rates about the X Y Z respectively
%             [quaternion_k_plus] = OMAS_axisTools.updateQuaternion(obj.VIRTUAL.quaternion,ENUrates,dt)
%             R_ENU = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);

            % [ TO BE CONFIRMED ] ENU ROTATIONS ARE CURRENTLY MAPPED TO THE
            % 'body rates' - globalRates 
            ENUrates = localENUState(omegaIndices); % rates about the X Y Z respectively
            globalRates = zeros(3,1);
            globalRates(1) = ENUrates(3);
            globalRates(2) = ENUrates(2);
            globalRates(3) = ENUrates(1);
            [quaternion_k_plus] = OMAS_axisTools.updateQuaternion(obj.VIRTUAL.quaternion,globalRates,dt);
            R_ENU = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = R_ENU*localENUState(velocityIndices);
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*globalVelocity_k_plus; % ASSUME INSTANTANEOUS
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.priorState = obj.localState;
            obj.localState = localENUState;                                % Reassign the obj.localstate
            % NO WAYPOINT CONSIDERATION
        end
        % UPDATE THE GLOBAL PROPERTIES ( LOCAL AXIS REMAINS STATIC )
        function [obj] = updateGlobalProperties_fixedFrame(obj,dt,localState_update)
           % This function computes the new global parameters for the given
           % agent based on its new state.
            
            % CONSTANT PARAMETERS
            velocityIndices = 7:9;
           
            % ////////// DEFINE PREVIOUS PARAMETERS ///////////////////////
            globalPosition_k = obj.VIRTUAL.globalPosition;                 % 3D although the 2D
            quaternion_k     = obj.VIRTUAL.quaternion;

            velocity_ENU_k_plus = localState_update(velocityIndices,1); 
                        
            % ////////// [WORKING: NON-ROTATED]
            quaternion_k_plus = quaternion_k;
            
            %NEW ROTATION MATRIX
            [R_update] = OMAS_axisTools.quaternionToRotationMatrix(quaternion_k_plus);
            localVelocityUpdate = velocity_ENU_k_plus;
            
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = R_update*localVelocityUpdate;
            globalPosition_k_plus = globalPosition_k + dt*globalVelocity_k_plus;
            
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.rotationMatrix = R_update;
            obj.localState = localState_update;  
        end
    end
    
    % ///////////////////// GENERAL STATIC METHODS ////////////////////////
    methods (Static)
        % ATTEMPT TO DEAL WITH NESTED CELL ARRAY INPUTS
        function [inputParameters] = inputHandler(inputParameters)
            % ATTEMPT TO REMOVE NESTED CELL ARRAYS
            if numel(inputParameters) == 1 && iscell(inputParameters)
                t = 1; % limit the number of tries..
                while numel(inputParameters) == 1 && t < 20
                    inputParameters = inputParameters{1};                  % Unpack nested input sets
                    t = t + 1;
                end
            end 
        end
        % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
        function [config] = configurationParser(defaultConfig,inputParameters)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called in the class constructor
                        
            assert(nargin == 2,'Parser must be provided with an object/structure and a cell parameter vector');
            
            % MOVE THROUGH THE PARAMETER PAIRS (i.e. ,'radius',1,...)
            for parameterIndex = 1:numel(inputParameters)
                % FOR EACH USER-DEFINED PARAMETER
                givenParameter = inputParameters{parameterIndex};
                if ~ischar(givenParameter)
                    continue
                else
                    % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                    if isstruct(defaultConfig) && isfield(defaultConfig,givenParameter)
                        defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};  % Make a substitution
                    end
                    if ~isstruct(defaultConfig) && isprop(defaultConfig,givenParameter)
                        defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};   % Make a substitution
                    end
                end
            end
            config = defaultConfig;
        end
        % BOUND VALUE
        function [boundedVector] = boundValue(inputVector,lower,upper)
            % This function simply bounds a value between two given values
            
            boundPerDim = size(lower,1) == size(inputVector,1);
            assert(size(lower,1) == size(upper,1),'Please specify the same number of upper and lower bounds');
                        
            if boundPerDim
                % ELSE ASSUME THERE IS UPPER AND LOWER BOUND FOR EACH VECTOR
                % DIMENSION
                for ind = 1:size(inputVector,1)
                    if inputVector(ind) < lower(ind)
                        inputVector(ind) = lower(ind);
                    end
                    if inputVector(ind) > upper(ind)
                        inputVector(ind) = upper(ind);
                    end
                end
            else
                % THERE IS ONE BOUND FOR ALL DIMENSIONS
                for ind = 1:size(inputVector,1)
                    if inputVector(ind) < lower
                        inputVector(ind) = lower;
                    end
                    if inputVector(ind) > upper
                        inputVector(ind) = upper;
                    end 
                end
            end
            boundedVector = inputVector;
        end
    end
    
    % //////////////////// SIMULATION STATIC METHODS //////////////////////
    methods (Static)
        % APPLY NAMING CONVENTION
        function [namestr]  = getObjectName(objectID)
            % No input name string specified, use greek naming scheme
            defaultID = { 'alpha', 'beta','gamma', 'delta','epsilon'  ...
                        ,  'zeta',  'eta','theta',  'iota','kappa'    ...
                        ,'lambda',   'mu',   'nu',    'xi','omicron'  ...
                        ,    'pi',  'rho','sigma',   'tau','upsilon'  ...
                        ,   'phi',  'chi',  'psi', 'omega'};
            % Generate default name string
            deflength = length(defaultID);
            cycle = ceil(objectID/deflength)-1;
            IDval = objectID - (cycle*deflength);                          %  start from 'alpha' and '001'
            namestr = strcat(defaultID{IDval},num2str(cycle+1,'%03d'));
        end
        % APPLY OBJECT-ID CONVENTION
        function [objectID] = getObjectID()
            % This function allocates a basic ID convention to the class
            % heirarchy to provide object automatic numeration.
            
            persistent objectCount;
            
            % DEFINE ID BASED ON EXISTING OBJECTS
            if isempty(objectCount)
                % INITIAL OBJECT
                objectCount = 1;
            else
                % OBJECTS EXIST
                objectCount = objectCount + 1;
            end
            % ALLOCATE OBJECT ID
            objectID = objectCount;
        end  
    end
end