%% GENERIC object CLASS (objectDefinition.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic object and import its 
% variables into the simulation space for the purpose of multi-vehicle 
% control simulation.

% Author: James A. Douthwaite 20/05/2016

classdef objectDefinition < handle
    %% OBJECT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % SIMULATION IDENTIFIERS
        objectID = 0;
        name;
        radius = 1;                 % Characteristic radius 
        % KINEMATIC PROPERTIES
        localState  = zeros(6,1);   % 6DOF state vector
        % GEOMETRIC PROPERTIES
        GEOMETRY = struct(...
            'vertices',[],...
            'faces',[],...
            'normals',[],...
            'centroid',zeros(3,1));             
    end
    % Private (OpenMAS-facing) properties
    properties (Access = private)
        % GLOBAL PROPERTIES (Virtual (SIMULATION) data container)
        GLOBAL; 
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods (Access = public)
        % Constructor
        function [this] = objectDefinition(varargin)
            % This function is designed to automatically apply a generic
            % naming convention to all objects generated prior to a
            % simulation.
            % INPUTS:
            % varargin - A cell array of property/value pairs 
            % OUTPUTS:
            % this     - The generated object
            
            % Check all system paths are available
            this.GetPathDependencies();
            
            % Assign the objectID (if required)
            if this.objectID == 0                                           % Generate objectID if not allocated by a super class
                this.objectID = this.CreateobjectID();                       % Assign the number as an ID tag
                spawnNewobject = 1;
            else
                spawnNewobject = 0;
            end
            
            % Assign a name string
            if numel(this.name) < 1   
                this.name = this.CreateName(this.objectID);                   % No input name string specified, use greek naming scheme
            end  
            
            % Assign default 
            this.GLOBAL   = this.CreateGLOBAL();    % The GLOBAL structure
            this.GEOMETRY = this.CreateGEOMETRY();  % Get the geometry of the object

            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
            
            % ////////// DISPLAY THE CURRENT object PROPERTIES ////////////
            if spawnNewobject
                fprintf('objectcount: %d\tname: %s\ttype: %s\n',this.objectID,this.name,class(this));
            end            
        end 
        % Setup (default 6-dof state)
        function [this] = setup(this,localXYZvelocity,localXYZrotations)   % [x y z phi theta psi]
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % INPUTS:
            
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi]

            % Standard state vector is the 6DOF: x_k = [x y z phi theta psi]
            this.localState = [zeros(3,1);localXYZrotations];
        end
        % Main  (default 6-dof state)
        function [this] = main(this,ENV,varargin)
            % This is a generic process cycle of an object that accepts no
            % input commands/feedback and simply updates its states based
            % on its current attributes
            % INPUT:
            % varargin - Generic variable input container, containing the simulation timestep
            % OUTPUT:
            % this - The updated object
            
            % DETERMINE THE TIMESTEP
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('object TIME packet is invalid.');
            end
            
            
            % The state differential
            x_dot = [0;0;0;0;0;0];
            
            % MAINTAIN LOCAL INPUT CONDITIONS
            % Default behaviour of a system is to move with constant input
            % conditions.
            [dXdt]   = this.SingleIntegratorDynamics(this.localState,x_dot);
            newState = this.localState + dt*dXdt;
            
            % UPDATE THE CLASS GLOBAL PROPERTIES
            this = this.GlobalUpdate(dt,newState);
        end
    end
    
    %% /////////////////////// CHECK FUNCTIONS ////////////////////////////
    methods
        % Get the 'idle' status
        function [flag] = IsIdle(this)
            flag = this.GetGLOBAL('idleStatus');
        end
        % Get the dimensionality of the object
        function [flag] = Is3D(this)
            flag = this.GetGLOBAL('is3D');
        end
    end
    %% ////////////////////// GET/SET FUNCTIONS ///////////////////////////
    methods
        % Set the radius
        function set.radius(this,radius)
            this.radius = radius;               % Set the radius
            this.SetGLOBAL('radius',radius);    % Set the simulation radius
            % If agent has GEOMETRY
            if ~isempty(this.GEOMETRY.vertices) 
                % Normalise the geometry
                this.GEOMETRY = OMAS_graphics.normalise(this.GEOMETRY);
                % Scale the geometry to the radius
                this.GEOMETRY = OMAS_graphics.scale(this.GEOMETRY,this.radius);
            end
        end
        % Get the GLOBAL structure
        function [param] = GetGLOBAL(this,label)
            % By default return the GLOBAL structure
            if nargin < 2
                param = this.GLOBAL;
                return
            end
            assert(ischar(label),"Expecting a string parameter label");
            % Return the associated parameter
            param = this.GLOBAL.(label);
        end  
        % Set the virtual structure
        function [this]  = SetGLOBAL(this,label,value)
            % By default set the GLOBAL structure
            if nargin < 3
                % Input sanity check
                assert(isstruct(label),"The object's GLOBAL field must be a structure.");
                % Set the GLOBAL structure
                this.GLOBAL = label;
                return
            end
            assert(ischar(label),"Expecting a string parameter label");
            % Assign the virtual structure
            this.GLOBAL.(label) = value;
        end        
    end
    %% ////////////////// BASIC PROPERTY ASSIGNMENTS //////////////////////
    methods
        % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
        function [this]     = ApplyUserOverrides(this,pairArray)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called in the class constructor
            
            % Input sanity check #1
            if nargin < 2 || numel(pairArray) == 0
                return
            end
            % Call the all parameter overrider
            [this] = GetParameterOverrides_recursive(this,pairArray);
        end
        % Get the object GEOMETRY structure
        function [GEOMETRY] = CreateGEOMETRY(this)
            % This function prepares object STL files for representation.
            % ASSUMPTIONS:
            % - The object stl file has the same name as the object being simulated.
            % - The STL is correctly rotated to match the axes of XYZ of the agents
            %   local frame.
            
            % This functions allow additional operations to the GEOMETRY
            % structure to happen other than sourcing the basic structure.
            
            % Call the generic object.stl geometry parser
            [GEOMETRY] = this.GetObjectGeometry(this);
        end
    end
    % Only the objectDefinition class has access
    methods (Static, Access = private)  
        % Get the object GLOBAL structure
        function [GLOBAL]   = CreateGLOBAL()
            % This is function that assembles the template desciption of an
            % object
            
            % Define the GLOBAL structure
            GLOBAL = struct();
            GLOBAL.type = OMAS_objectType.misc;
            GLOBAL.hitBoxType = OMAS_hitBoxType.none;
            GLOBAL.radius = 0.5;                        % Diameter of 1m
            GLOBAL.colour = rand(1,3,'single');         % Virtual colour (for plotting)
            GLOBAL.symbol = 'square';                   % Representative symbol
            GLOBAL.position = [0;0;0];                  % Global Cartesian position
            GLOBAL.velocity = [0;0;0];                  % Global Cartesian velocity
            GLOBAL.quaternion = [1;0;0;0];              % Global quaternion pose
            GLOBAL.idleStatus = true;                   % object idle logical
            GLOBAL.is3D = true;                         % object operates in 3D logical
            GLOBAL.priorState = [];                            
        end
        % Get the object name
        function [namestr]  = CreateName(objectID)
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
        % Get the object ID
        function [objectID] = CreateobjectID()
            % This function allocates a basic ID convention to the class
            % heirarchy to provide object automatic numeration.
            
            persistent objectCount;
            
            % DEFINE ID BASED ON EXISTING objectS
            if isempty(objectCount)
                % INITIAL object
                objectCount = 1;
            else
                % objectS EXIST
                objectCount = objectCount + 1;
            end
            % ALLOCATE object ID
            objectID = uint16(objectCount);
        end 
    end

    %% STATE UPDATE METHODS
    methods (Static)
        % STATE UPDATE - TWO WHEELED DIFFERENTIAL ROBOT
        function [dXdt] = DifferencialDriveDynamics(X,d,a_left,a_right)
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
        function [dXdt] = UnicycleDynamics(X,speed,headingRate)
            % This function provides the dynamics of a uni-cycle model
            % moving in on a 2D plane.
                        
            % THE STATE UPDATE
            dXdt = zeros(3,1);
            dXdt(1) = speed*cos(X(3));
            dXdt(2) = speed*sin(X(3));
            dXdt(3) = headingRate;
        end
        % STATE UPDATE - DOUBLE INTEGRATOR
        function [dXdt] = DoubleIntegratorDynamics(X,U)
            % The relationship between the input and state vector is by
            % definition of a double integrator: dX(dX/dt)/dt = U
            
            % X(t) here is assumed to be of the form X(t) = [q;dq]
            
            qDim = length(X)/2;
            
            % dXdt is of the form dX(t) = [dq;ddq];
            dXdt = [X(qDim+1:2*qDim,1);     % Velocity states
                    U];                     % Acceleration states
        end
        % STATE UPDATE - SINGLE INTEGRATOR
        function [dXdt] = SingleIntegratorDynamics(X,U)
            % The relationship between the input and state vector is by
            % definition of a single integrator: dX/dt = U
            
            % Sanity check
            assert(length(X) == length(U),"Expecting the number of inputs to correspond to the state vector.");
            
            
            
            % X(t) here is assumed to be of the form X(t) = [q]
            
            % The change in state is literally the input
            dXdt = U;            
        end
        % STATE UPDATE - SIMPLE EULER VELOCITIES
        function [dXdt] = SimpleDynamics(X,velocity_k_plus,omega_k_plus)
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
    end
    %% GLOBAL REPRESENTATION UPDATE METHODS
    methods
        % INITIAL 3D STATE - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
        function [this] = setup_3DVelocities(this,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y z phi theta psi; dx dy dz dphi dtheta dpsi]
            
            % BUILD STATE VECTOR
            this.localState = zeros(12,1);
            this.localState(4:5) = localXYZrotations(1:2);
            %this.localState(4:6) = localXYZrotations;       % Represent yaw relative to initial position    
            this.localState(7:9) = localXYZVelocity;
            % Define the 
            this.SetGLOBAL('priorState',this.localState);
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x_t;x_dot]'
        function [this] = GlobalUpdate_3DVelocities(this,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y z phi theta psi dx dy dz dphi dtheta dpsi].
            
            positionIndices = 1:3;
            angleIndices = 4:6;
            
            % USE THE 'RATE' STATES DIRECTLY
            localLinearRates  = eulerState(positionIndices+6);
            localAngularRates = eulerState(angleIndices+6);
            % Get current properties
            p_k = this.GetGLOBAL('position'); 
            v_k = this.GetGLOBAL('velocity'); 
            q_k = this.GetGLOBAL('quaternion');   
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION % [ IF THE QUATERNION DEFINES GB ]

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
            this = this.GlobalUpdate_direct(...
                p_k_plus,...    % Global position at k plius
                v_k_plus,...    % Global velocity at k plus
                q_k_plus);      % Quaternion at k plus
        end
        % GLOBAL UPDATE - EULER 6DOF(3DOF) [x y z phi theta psi],[x y psi]
        function [this] = GlobalUpdate(this,dt,eulerState)
            % This function computes the position, velocity and 
            % quaternion in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % This function is intended for state vectors of the form:
            % - [x y z phi theta psi]
            % - [x y phi]
                                   
            % NOTATION INDICIES
            if this.Is3D
                positionIndices = 1:3;
                eulerIndices = 4:6;     % The rotations about X,Y,Z
            else
                positionIndices = 1:2;
                eulerIndices = 3;       % The rotations about Z
            end
          
            % Get current properties
            p_k = this.GetGLOBAL('position'); 
            v_k = this.GetGLOBAL('velocity'); 
            q_k = this.GetGLOBAL('quaternion');  
            x_k = this.GetGLOBAL('priorState');  
            
            % Populate the prior state if necessary
            if isempty(x_k)
                x_k = this.localState;
            end  
            
            % Equivalent rates
            localLinearRates  = (eulerState(positionIndices) - x_k(positionIndices))/dt;
            localAngularRates = (eulerState(eulerIndices) - x_k(eulerIndices))/dt;  
            
            % EULER ROTATIONS -> GLOBAL ROTATIONS
            if numel(localAngularRates) ~= 3
                localAngularRates = [0;0;1]*localAngularRates; 
                localLinearRates  = [localLinearRates;0];
            end
                        
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(q_k');
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*localAngularRates;
            % UPDATE THE QUATERNION POSE
            q_k_plus = OMAS_geometry.integrateQuaternion(q_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(q_k_plus');    
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            velocity_k_plus = R_k_plus'*localLinearRates;
            position_k_plus = p_k + dt*v_k;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            this = this.GlobalUpdate_direct(...
                position_k_plus,... 	% Global position at k plius
                velocity_k_plus,...     % Global velocity at k plus
                q_k_plus);              % Quaternion at k plus
        end
        % GLOBAL UPDATE - EULER 6DOF(3DOF) (NO ROTATIONS)
        function [this] = GlobalUpdate_fixedFrame(this,dt,eulerState)
            % This function computes the position and velocity whilst 
            % maintaining a fixed quaternion rotation (constant reference 
            % orientation) in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % NOTATION INDICIES
            if this.Is3D
                positionIndices = 1:3;
            else
                positionIndices = 1:2;
            end
            
            % Get current properties
            position_k   = this.GetGLOBAL('position'); 
            velocity_k   = this.GetGLOBAL('velocity'); 
            quaternion_k = this.GetGLOBAL('quaternion');  
            state_k      = this.GetGLOBAL('priorState');
            
            % Equivalent velocity
            velocity_k_plus = (eulerState(positionIndices) - state_k(positionIndices))/dt;
            
            % ROTATION RATES ABOUT THE GLOBAL AXES
            if numel(eulerState) ~= 6
                velocity_k_plus = [velocity_k_plus;0];
            end
            
            % NEW ROTATION MATRIX (G>B)            
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(quaternion_k);
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            velocity_k_plus = R_k_plus'*velocity_k_plus;
            position_k_plus = position_k + dt*velocity_k;
            
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.GlobalUpdate_direct(...
                position_k_plus,...     % Global position at k plius
                velocity_k_plus,...     % Global velocity at k plus
                quaternion_k_plus,...   % Quaternion at k plus
                eulerState);            % The new state for reference
        end
        % GLOBAL UPDATE - DIRECTLY (global default)
        function [this] = GlobalUpdate_direct(this,p,v,q)
            % Under this notation, the state vector already contains the
            % global parameters of the object.
            % INPUTS:
            % globalPosition - 3D global cartesian position given in the ENU coordinate frame.
            % globalVelocity - 3D global cartesian velocity given in the ENU coordinate frame. 
            % quaternion     - The new quaternion pose of body in the global frame.
            % R              - The rotation of the body
            % this.localState - The previous localState (independant of convention)
            % eulerState     - The new state as reported by the agent            
            
            % Input sanity check
            assert(IsColumn(p,3),'Global position must be a 3D column vector [3x1].');
            assert(IsColumn(v,3),'Global velocity must be a 3D column vector [3x1].');
            assert(IsColumn(q,4),'Global pose must be a 4D quaternion vector [4x1].');
            assert(size(this.localState,2) == 1,'The length of the objects state update must match the its local state.');
            
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            % Convert the quaternion to the equivalent rotation matrix
            R = OMAS_geometry.quaternionToRotationMatrix(q);
            % Assign the global parameters
            this.SetGLOBAL('position',p);                 	% Reassign the global position
            this.SetGLOBAL('velocity',v);                 	% Reassign the global velocity
            this.SetGLOBAL('quaternion',q);                	% Reassign the quaternion
            this.SetGLOBAL('R',R);                          % New rotation matrix
            this.SetGLOBAL('priorState',this.localState);  	% Record the previous state
        end
    end
    methods (Static)
        % Generic .STL to geometry parser
        function [GEOMETRY] = GetObjectGeometry(entity)
            % This function prepares object STL files for representation.
            % ASSUMPTIONS:
            % - The object stl file has the same name as the object being simulated.
            % - The STL is correctly rotated to match the axes of XYZ of the agents
            %   local frame.
            
            % ABSOLUTE PATH TO THE object STL LOCATION
            absPath = mfilename('fullpath');
            absPath = strrep(absPath,'objectDefinition','');               % Get the current m-files directory
            
            % GET THE AGENT SUPERCLASSES
            aliases = superclasses(entity);
            aliases = vertcat(class(entity),aliases);
            % MOVE THROUGH THE CLASS HEIRARCHY
            for aliasNo = 1:numel(aliases)
                % BUILD THE STL ASSOCIATION
                filename = strcat(absPath,char(aliases(aliasNo)),'.stl');  % Build the associated STL path
                % GET THE object PATCH FROM STL FILE (in the object directory)
                [GEOMETRY,getFlag] = OMAS_graphics.importStlFromFile(filename);
                if getFlag
                    break                                                  % Successful import
                end
            end
            
            % IF PATCH RETURNED, PROCESS FOR SIMULATION
            if isstruct(GEOMETRY)
                % Get the entities representative radius
                entityRadius = entity.GetGLOBAL('radius');
                % NORMALISE THE IMPORTED GEOMETRY
                GEOMETRY = OMAS_graphics.normalise(GEOMETRY);              % Normalise
                GEOMETRY = OMAS_graphics.scale(GEOMETRY,entityRadius);     % Scale
                GEOMETRY.normals = OMAS_graphics.normals(GEOMETRY);
                % ADD CENTROID
                GEOMETRY.centroid = zeros(1,3);                            % Assume vertices are relative to a centroid
            else
                % GEOMETRIC PROPERTIES
                GEOMETRY = struct('vertices',[],...
                                  'faces',[],...
                                  'normals',[],...
                                  'centroid',zeros(1,3));                  % If the object is something other than a point.        
            end
            % Return patch, empty if not successful
        end
        % Dependancy check for all objects
        function GetPathDependencies()
            % This program preforms a check of the dependencies necessary
            % for the instantiation of object classes.
            
            % Get the path to the install directory
            repoString = mfilename('fullpath');
            repoString = repoString(1:(strfind(repoString,'objects')-1));
            % Add the path to the repo's environment 
            addpath([repoString,'environment']);
            % Check all the OMAS_dependencies
            if ~OMAS_system.GetFileDependancies()
                error('Unable to get OpenMAS dependancies.');
            end
        end
    end
end