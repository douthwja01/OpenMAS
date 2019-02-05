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
        localState  = zeros(12,1);                                         % Other properties derived from the state.
        % GEOMETRIC PROPERTIES
        GEOMETRY = struct('vertices',[],...
                          'faces',[],...
                          'normals',[],...
                          'centroid',zeros(3,1));                          % If the object is something other than a point.
        % VIRTUAL PROPERTIES (Virtual (SIMULATION) data container)
        VIRTUAL = struct();                                                % Object operates in 3D logical
    end
    % CLASS METHODS
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
            
            % ///////////////////// OBJECT SETUP //////////////////////////
            % ASSIGN VIRTUAL PROPERTIES
            [obj.VIRTUAL] = obj.getObjectVIRTUAL();
            % GET GEOMETRY IF POSSIBLE (FOR ALL CLASS SIBLINGS)
            [obj.GEOMETRY] = obj.getObjectGeometry(obj);                   % Get the geometry of the object
            % ALLOCATE OBJECT ID
            spawnNewObject = 0;
            if isempty(obj.objectID)                                       % Generate objectID if not allocated by a super class
                spawnNewObject = 1;
                [obj.objectID] = obj.getObjectID();                        % Assign the number as an ID tag
            end
            % GENERATE OBJECT NAME IF REQUIRED
            if isempty(obj.name) || numel(obj.name) < 1   
                obj.name = obj.getObjectName(obj.objectID);                % No input name string specified, use greek naming scheme
            end            
            
            % ////////// HANDLING THE USERS INPUTS AS OVERRIDES ///////////
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
            
            % ////////// DISPLAY THE CURRENT OBJECT PROPERTIES ////////////
            if spawnNewObject
                fprintf('Objectcount: %d\tname: %s\ttype: %s\n',obj.objectID,obj.name,class(obj));
            end            
        end 
        % ///////////////////// SETUP FUNCTION ////////////////////////////
        function [obj] = setup(obj,localXYZvelocity,localXYZrotations)     % [x y z phi theta psi]
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % INPUTS:
            
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi]
            
            % BUILD STATE VECTOR
            obj.localState = zeros(6,1);
            obj.localState(4:6) = localXYZrotations;
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj.VIRTUAL.priorState = obj.localState;
        end
        % ////////////////// THE DEFAULT OBJECT CYCLE /////////////////////
        function [obj] = main(obj,TIME,varargin)
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
            % Default behaviour of a system is to move with constant input
            % conditions.
            [dXdt] = obj.dynamics_singleIntegrator(obj.localState,[0;0;0],[0;0;0]);
            newState = obj.localState + dt*dXdt;
            
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties_ENU(dt,newState);
        end
    end
    %% ///////////////// BASIC STATE UPDATE FUNCTIONS /////////////////////
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
        
        
        % UPDATE AGENT STATE VECTOR FROM LOCAL ACCELERATIONS
        function [dXdt] = dynamics_doubleIntegrator(obj,X,linearAcceleration,headingAcceleration)
            % This function provides a basic model for the agent
            % kinematics, controlled by linear and angular acceleration inputs            
            
            % ASSUMPTIONS
            % - The state vector is of the form:
            %   X = [x eta ; dx deta ] ; U = [ddx ; ddeta]
            % - The inputs are an order higher than the highest order
            %   states.
            
            % SET THE INPUT ORDER
            % The input "linearAcceleration" and "headingAcceleration" are
            % related to the second order of the states directly.
            systemOrder = 2;
            % PASS THE SECOND ORDER STATES TO THE N-ORDER INTEGRATOR
            dXdt = obj.dynamics_NIntegrator(X,systemOrder,linearAcceleration,headingAcceleration);
        end
        % UPDATE AGENT STATE VECTOR FROM LOCAL VELOCITIES
        function [dXdt] = dynamics_singleIntegrator(obj,X,linearRates,headingRates)
            % This function preforms the state update directly from a
            % velocity vector in the body axis system.
            
            % ASSUMPTIONS
            % - The state vector is of the form:
            %   X = [x eta] ; U = [dx ; deta]
            % - The inputs are an order higher than the highest order
            %   states.
            
            % SET THE INPUT ORDER
            % The input "linearRates" and "headingRate" are related are the
            % first order of the states directly.
            systemOrder = 1; 
            % PASS THE FIRST ORDER STATES TO THE N-ORDER INTEGRATOR
            dXdt = obj.dynamics_NIntegrator(X,systemOrder,linearRates,headingRates);
        end
    end
    
    methods (Static)
%         % DYNAMICS - DOUBLE INTEGRATOR
%         function [dXdt] = dynamics_doubleIntegrator(X,a_linear,a_angular)
%             
%             if numel(a_angular) == 1
%                 is3D = 0;
%             else
%                 is3D = 1;
%             end
%             
%             numStates = numel(X);           % The number of states
%             stateRepetitions = 2;
%             
%             assert(mod(numel(X),stateNum) == 0,'The number of states does not correspond to the input order.');
%              
%             length_q = numStates/stateRepetitions;
%             
%             U = [a_linear;a_angular];       % Concatinate inputs
%             
%             
%             zeroMatrix = zeros(numStates);
%             
%             dXdt = [zeroMatrix,       eye(numStates);
%                     zeroMatrix,zeroMatrix(numStates)]*X + [0;1]*U;   % Double integrator dynamics
%         end
        % STATE UPDATE - TWO WHEELED DIFFERENTIAL ROBOT
        function [dXdt] = dynamics_differencialDrive(X,d,a_left,a_right)
            % Becker, M. (2006). Obstacle avoidance procedure for mobile robots. ABCM Symposium Series
            % [x;y;theta;v;omega]

            % BUILD THE STATE DIFFERENTAL
            dXdt = zeros(5,1);
            dXdt(1) = X(4)*cos(X(3));             % x velocity
            dXdt(2) = X(4)*sin(X(3));             % y velocity
            dXdt(3) = X(5);                       % rotational rate
            dXdt(4) = (a_right + a_left)/2;       % Linear acceleration
            dXdt(5) = (a_right - a_left)/d;       % Angular acceleration
        end
        % STATE UPDATE - UNI-CYCLE ROBOT
        function [dXdt] = dynamics_unicycle(X,speed,headingRate)
            % This function provides the dynamics of a uni-cycle model
            % moving in on a 2D plane.
                        
            % THE STATE UPDATE
            dXdt = zeros(3,1);
            dXdt(1) = speed*cos(X(3));
            dXdt(2) = speed*sin(X(3));
            dXdt(3) = headingRate;
        end
        % STATE UPDATE - SIMPLE EULER VELOCITIES
        function [dXdt] = dynamics_simple(X,velocity_k_plus,omega_k_plus)
            % This function assumes that the velocity changes are
            % implemented this timestep directly, integration then occurs
            % including the updates from this timestep:
            % State vector is assumed:
            % 2D - [x y psi dx dy dpsi]'
            % 3D - [x y z phi theta psi dx dy dz dphi dtheta dpsi]'
            
            % CHECK EXPECTED INPUT DIMENSIONS
            Aflag = numel(velocity_k_plus) == 3 && numel(omega_k_plus) == 3; % IS 3D
            Bflag = numel(velocity_k_plus) == 2 && numel(omega_k_plus) == 1; % IS 2D
            assert((Aflag || Bflag) == 1,'Integrator inputs are of the wrong dimensions');
            
            % Either rotation occurs in one dimension (2D) or in three (3D)
            rotationDOF = numel(omega_k_plus); % Is heading a yaw (2D) or [phi,theta psi] (3D)
            if rotationDOF == 1
                linearIndices = 1:2;
                angularIndices = 3;
            else
                linearIndices = 1:3;
                angularIndices = 4:6;
            end
            
            % THE STATE DIFFERENTIAL
            dXdt = zeros(numel([linearIndices,angularIndices]),1);
            dXdt(linearIndices) = velocity_k_plus;
            dXdt(angularIndices) = omega_k_plus;
        end
        % STATE UPDATE - N-ORDER INTEGRATOR; INPUT OF THE HIGHEST ORDER
        function [dXdt] = dynamics_NIntegrator(X,order,U_linear,U_angular)
            % This function preforms the state update from an input vector
            % of relation N to the given state parameter.
            % ASSUMPTIONS:
            % - The agent is capable of changing its heading instanaeously.
            % - obj.localState is of length N*states of the form [x,dx,ddx,...]
            
            % INPUTS:
            % X            - The current state of order 0:n-1
            % n            - The order of the input
            % linearRates  - Defined as the rates along the body axes (m/s)
            % headingRates - Defined as the euler rates [dphi;dtheta;dpsi] (rad/s)
            
            % CHECK EXPECTED INPUT DIMENSIONS
            Aflag = numel(U_linear) == 3 && numel(U_angular) == 3;   % IS 3D
            Bflag = numel(U_linear) == 2 && numel(U_angular) == 1;   % IS 2D
            assert((Aflag + Bflag) == 1,'Integrator inputs are of the wrong dimensions');
            
            % CALCULATE THE STATE-INPUT PARAMETERS
            % Either rotation occurs in one dimension (2D) or in three (3D)
            rotationDOF = numel(U_angular); % Is heading a yaw (2D) or [phi,theta psi] (3D)
            if rotationDOF == 1
                linearIndices = 1:2;
                angularIndices = 3;     % Is 2D
            else
                linearIndices = 1:3;
                angularIndices = 4:6;   % Is 3D
            end
            
            % 'n' represents the order of the inputs, and also an order beyond 
            % the state vector (i.e X = [x dx], U = [ddx])
            
            % ENSURE THIS CONVENTION
            stateNum = numel([linearIndices,angularIndices]);
            assert(mod(numel(X),stateNum) == 0,'The number of states does not correspond to the input order.');
             
            % THE 'n+1'th ORDER INPUTS
            U = zeros(numel(X),1);
            U(linearIndices,1)  = U_linear;
            U(angularIndices,1) = U_angular;                     % The 'nth' order axial inputs
                        
            % MAP THE HIGHER ORDER STATES TO THE LOWER STATE CHANGES
            % dX = [0 1 0][x;dx]
            %      [0 0 1][U]
            % CONSTRUCT MAPPING MATRICES
            mappingMatrix = [zeros(stateNum*(order),stateNum),diag(ones(stateNum*(order),1))];
            
            
            assert(size(mappingMatrix,2) == size([X;U],1),'Integration dimensions are incorrect');
            % INTEGRATE THE LOWER ORDER STATES, CONCATINATE THE INPUTS AS
            dXdt = mappingMatrix*[X;U];  
        end
    end
    
    %% ////////////////// SIMULATION & CORE INTERFACES ////////////////////
    methods
        % INITIAL 3D STATE - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
        function [obj] = initialise_3DVelocities(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            
            % BUILD STATE VECTOR
            obj.localState = zeros(12,1);
            obj.localState(4:6) = localXYZrotations;
            obj.localState(7:9) = localXYZVelocity;
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj.VIRTUAL.priorState = obj.localState;
        end
        
        % GLOBAL UPDATE - DEFAULT 
        function [obj] = updateGlobalProperties(obj,dt,eulerState)
            % SIMPLY A MAP TO THE 3D UPDATE FUNCTION
            obj = obj.updateGlobalProperties_ENU(dt,eulerState);
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x_t;x_dot]'
        function [obj] = updateGlobalProperties_3DVelocities(obj,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y z phi theta psi dx dy dz dphi dtheta dpsi].
            
            positionIndices = 1:3;
            angleIndices = 4:6;
            
            % USE THE 'RATE' STATES DIRECTLY
            velocity_k_plus = eulerState(positionIndices+6);
            localAxisRates  = eulerState(angleIndices+6);
            
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
        % GLOBAL UPDATE - EULER 6DOF(3DOF) [x y z phi theta psi],[x y psi]
        function [obj] = updateGlobalProperties_ENU(obj,dt,eulerState)
            % This function computes the position, velocity and 
            % quaternion in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % This function is intended for state vectors of the form:
            % - [x y z phi theta psi]
            % - [x y phi]
                                   
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
                eulerIndices = 4:6;                     % The rotations about X,Y,Z
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
                eulerIndices = 3;                       % The rotations about Z
            else
                error('State notation not recognised');
            end
            
            % EQUIVALENT RATES
            velocity_k_plus   = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt;
            eulerRates_k_plus = (eulerState(eulerIndices) - obj.VIRTUAL.priorState(eulerIndices))/dt;  
            
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isprop(obj,'targetWaypoint')
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj.VIRTUAL.idleStatus = 1;
                    velocity_k_plus   = zeros(numel(positionIndices),1);   % Freeze the agent
                    eulerRates_k_plus = zeros(numel(eulerIndices),1);      % Freeze the agent
                end   
            end

            % EULER ROTATIONS -> GLOBAL ROTATIONS
            if numel(eulerRates_k_plus) ~= 3
                localAxisRates = [0;0;1]*eulerRates_k_plus; 
                velocity_k_plus = [velocity_k_plus;0];
            else
                % CALCULATE THE GLOBAL AXIS RATES (FROM 3DOF)
                localAxisRates = eulerRates_k_plus;                        % Rotations rates are negative
            end
                        
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                % [ IF THE QUATERNION DEFINES GB ]
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
        % GLOBAL UPDATE - EULER 6DOF(3DOF) (NO ROTATIONS)
        function [obj] = updateGlobalProperties_ENU_fixedFrame(obj,dt,eulerState)
            % This function computes the position and velocity whilst 
            % maintaining a fixed quaternion rotation (constant reference 
            % orientation) in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
            else
                error('State notation not recognised');
            end
            
            % EQUIVALENT RATES
            velocity_k_plus = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt;
            
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isprop(obj,'targetWaypoint')
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj.VIRTUAL.idleStatus = 1;
                    velocity_k_plus = zeros(numel(positionIndices),1);     % Freeze the agent
                end   
            end
            
            % ROTATION RATES ABOUT THE GLOBAL AXES
            if numel(eulerState) ~= 6
                velocity_k_plus = [velocity_k_plus;0];
            end
            
            % NEW ROTATION MATRIX (G>B)            
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(obj.VIRTUAL.quaternion);
         
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end
        % GLOBAL UPDATE - DIRECTLY
        function [obj] = updateGlobalProperties_direct(obj,p,v,q,X)
            % Under this notation, the state vector already contains the
            % global parameters of the object.
            % INPUTS:
            % globalPosition - 3D global cartesian position given in the ENU coordinate frame.
            % globalVelocity - 3D global cartesian velocity given in the ENU coordinate frame. 
            % quaternion     - The new quaternion pose of body in the global frame.
            % R              - The rotation of the body
            % obj.localState - The previous localState (independant of convention)
            % eulerState     - The new state as reported by the agent            
            
            % SANITY CHECKS
            assert(numel(p) == 3 && size(p,2) == 1,'Global position must be a 3D column vector [3x1].');
            assert(numel(v) == 3 && size(v,2) == 1,'Global velocity must be a 3D column vector [3x1].');
            assert(numel(q) == 4 && size(q,2) == 1,'Global pose must be a 4D quaternion vector [4x1].');
%             assert(numel(X) == numel(obj.localState) && size(obj.localState,2) == 1,'The length of the objects state update must match the its local state.');
            
            % ///////////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalPosition = p;                                % Reassign the global position
            obj.VIRTUAL.globalVelocity = v;                                % Reassign the global velocity
            obj.VIRTUAL.quaternion = q;                                    % Reassign the quaternion
            obj.VIRTUAL.R = OMAS_geometry.quaternionToRotationMatrix(q);   % K_plus rotation
            obj.VIRTUAL.priorState = obj.localState;                       % Record the previous state
            obj.localState = X;                                            % Update the current state
        end
    end
    % //////////////////// SIMULATION STATIC METHODS //////////////////////
    methods (Static)
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
        % APPLY OMAS CONVENTION
        function [VIRTUAL]  = getObjectVIRTUAL()
            % This is function that assembles the template desciption of an
            % object
            % DEFINE THE VIRTUAL STRUCTURE
            VIRTUAL = struct(...
                         'type',OMAS_objectType.misc,...
                         'hitBoxType',OMAS_hitBoxType.none,...
                         'radius',0.5,...                                  % Diameter of 1m
                         'colour',rand(1,3,'single'),...                   % Virtual colour (for plotting)
                         'symbol','square',...                             % Representative symbol
                         'globalPosition',[0;0;0],...                      % Global Cartesian position
                         'globalVelocity',[0;0;0],...                      % Global Cartesian velocity
                         'quaternion',[1;0;0;0],...                        % Global quaternion pose
                         'idleStatus',logical(true),...                    % Object idle logical
                         'is3D',logical(true));                            % Object operates in 3D logical
        end
        % IMPORT THE CLASS GEOMETRY FILE
        function [geometry] = getObjectGeometry(obj)
            % This function prepares object STL files for representation.
            % ASSUMPTIONS:
            % - The object stl file has the same name as the object being simulated.
            % - The STL is correctly rotated to match the axes of XYZ of the agents
            %   local frame.
            
            % ABSOLUTE PATH TO THE OBJECT STL LOCATION
            absPath = mfilename('fullpath');
            absPath = strrep(absPath,'objectDefinition','');               % Get the current m-files directory
            
            % GET THE AGENT SUPERCLASSES
            aliases = superclasses(obj);
            aliases = vertcat(class(obj),aliases);
            % MOVE THROUGH THE CLASS HEIRARCHY
            for aliasNo = 1:numel(aliases)
                % BUILD THE STL ASSOCIATION
                filename = strcat(absPath,char(aliases(aliasNo)),'.stl');  % Build the associated STL path
                % GET THE OBJECT PATCH FROM STL FILE (in the object directory)
                [geometry,getFlag] = OMAS_graphics.importStlFromFile(filename);
                if getFlag
                    break                                                  % Successful import
                end
            end
            
            % IF PATCH RETURNED, PROCESS FOR SIMULATION
            if isstruct(geometry)
               % NORMALISE THE IMPORTED GEOMETRY
                [geometry] = OMAS_graphics.normalise(geometry);                 % Normalise
                [geometry] = OMAS_graphics.scale(geometry,obj.VIRTUAL.radius);  % Scale
                % ADD CENTROID
                geometry.centroid = zeros(1,3);                                 % Assume vertices are relative to a centroid
            else
                % GEOMETRIC PROPERTIES
                geometry = struct('vertices',[],...
                                  'faces',[],...
                                  'normals',[],...
                                  'centroid',zeros(3,1));                  % If the object is something other than a point.        
            end
            % Return patch, empty if not successful
        end
        % APPLY NAMING CONVENTION
        function [namestr]  = getObjectName(objectID)
            % No input name string specified, use greek naming scheme
            defaultID = { 'alpha', 'beta','gamma', 'delta','epsilon'  ...
                        ,  'zeta',  'eta','theta',  'iota','kappa'    ...
                        ,'lambda',   'mu',   'nu',    'xi','omicron'  ...
                        ,    'pi',  'rho','sigma',   'tau','upsilon'  ...
                        ,   'phi',  'chi',  'psi', 'omega'}; 
            % GENERATE DEFAULT NAME STRING
            defLength = length(defaultID);
            cycle = floor(double(objectID)/defLength);
            IDvalue = mod(objectID,defLength);
            if IDvalue == 0
                IDvalue = 24;
            end
            namestr = sprintf('%s%03d',defaultID{IDvalue},cycle+1);
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
            objectID = uint8(objectCount);
        end  
    end
end