%% GENERIC object CLASS (objectDefinition.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic object and import its variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite 20/05/2016

classdef objectDefinition < handle
    % object BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % SIMULATION IDENTIFIERS
        objectID = 0;
        name;
        % KINEMATIC PROPERTIES
        localState  = zeros(6,1);                                          % Other properties derived from the state.
        % GEOMETRIC PROPERTIES
        GEOMETRY = struct('vertices',[],...
                          'faces',[],...
                          'normals',[],...
                          'centroid',zeros(3,1));                          % If the object is something other than a point.
    %end
    %properties (Access = private)
        % VIRTUAL PROPERTIES (Virtual (SIMULATION) data container)
        VIRTUAL;                                  % object operates in 3D logical
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
            this.VIRTUAL  = this.CreateVIRTUAL();   % The VIRTUAL structure
            this.GEOMETRY = this.CreateGEOMETRY();  % Get the geometry of the object

            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
            
            % ////////// DISPLAY THE CURRENT object PROPERTIES ////////////
            if spawnNewobject
                fprintf('objectcount: %d\tname: %s\ttype: %s\n',this.objectID,this.name,class(this));
            end            
        end 
        % Setup (global default)
        function [this] = setup(this,localXYZvelocity,localXYZrotations)     % [x y z phi theta psi]
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % INPUTS:
            
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi]
            
            % Build state vector
            this.localState = zeros(6,1);
            this.localState(4:6) = localXYZrotations;
        end
        % Main (global default)
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
            
            % MAINTAIN LOCAL INPUT CONDITIONS
            % Default behaviour of a system is to move with constant input
            % conditions.
            [dXdt]   = this.dynamics_singleIntegrator(this.localState,[0;0;0;0;0;0]);
            newState = this.localState + dt*dXdt;
            
            % UPDATE THE CLASS GLOBAL PROPERTIES
            this = this.updateGlobalProperties(dt,newState);
        end
    end
    
    %% /////////////////////// CHECK FUNCTIONS ////////////////////////////
    methods
        % Get the 'idle' status
        function [flag] = IsIdle(this)
            flag = this.GetVIRTUALparameter('idleStatus');
        end
        % Get the dimensionality of the object
        function [flag] = Is3D(this)
            flag = this.GetVIRTUALparameter('is3D');
        end
    end
    %% ////////////////////// GET/SET FUNCTIONS ///////////////////////////
    methods
        % Get the VIRTUAL structure
        function [VIRTUAL]  = GetVIRTUAL(this)
            VIRTUAL = this.VIRTUAL;
        end  
        % Set the virtual structure
        function [this]     = SetVIRTUAL(this,VIRTUAL)
            % Input sanity check
            assert(isstruct(VIRTUAL),"The object's VIRTUAL field must be a structure.");
            % Assign the virtual structure
            this.VIRTUAL = VIRTUAL;
        end        
        % Get the virtual parameter
        function [value]    = GetVIRTUALparameter(this,label)
            % Input sanity check
            assert(ischar(label),"The property must be specified as a string.");
            assert(isprop(this,'VIRTUAL'),"object has no VIRTUAL property."); 
            assert(isfield(this.VIRTUAL,label),sprintf('%s is not a VIRTUAL parameter.',label));
            % Set the parameter to the field
            value = this.VIRTUAL.(label);
        end
        % Set the virtual parameter
        function [this]     = SetVIRTUALparameter(this,label,value)
            % Input sanity check
            assert(ischar(label),"The property must be specified as a string.");
            assert(isprop(this,'VIRTUAL'),"object has no VIRTUAL property."); 
            assert(isfield(this.VIRTUAL,label),sprintf('"%s" is not a VIRTUAL parameter.',label));
            % Set the parameter to the field
            this.VIRTUAL.(label) = value;
        end
        % Get general parameter
        function [value]	= GetParameter(this,label)
            % Input sanity check
            assert(ischar(label) && isprop(this,label),'The provided must be a parameter.');
            value = this.(label);
        end
        % Set general parameter
        function [this]     = SetParameter(this,label,value)
            % Input sanity check
            assert(ischar(label) && isprop(this,label),'The provided must be a parameter.');
            this.(label) = value;
        end
    end
    %% ////////////////// BASIC PROPERTY ASSIGNMENTS //////////////////////
    methods
        % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
        function [this]      = ApplyUserOverrides(this,pairArray)
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
        % Get the object VIRTUAL structure
        function [VIRTUAL]  = CreateVIRTUAL()
            % This is function that assembles the template desciption of an
            % object
            
            % Define the VIRTUAL structure
            VIRTUAL = struct();
            VIRTUAL.type = OMAS_objectType.misc;
            VIRTUAL.hitBoxType = OMAS_hitBoxType.none;
            VIRTUAL.radius = 0.5;                      % Diameter of 1m
            VIRTUAL.colour = rand(1,3,'single');       % Virtual colour (for plotting)
            VIRTUAL.symbol = 'square';                 % Representative symbol
            VIRTUAL.globalPosition = [0;0;0];          % Global Cartesian position
            VIRTUAL.globalVelocity = [0;0;0];          % Global Cartesian velocity
            VIRTUAL.quaternion = [1;0;0;0];            % Global quaternion pose
            VIRTUAL.idleStatus = true;                 % object idle logical
            VIRTUAL.is3D = true;                       % object operates in 3D logical
            VIRTUAL.priorState = [];                            
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

    %% ///////////////// BASIC STATE UPDATE FUNCTIONS /////////////////////
    methods (Static)
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
        % STATE UPDATE - DOUBLE INTEGRATOR
        function [dXdt] = dynamics_doubleIntegrator(X,U)
            % The relationship between the input and state vector is by
            % definition of a double integrator: dX(dX/dt)/dt = U
            
            % X(t) here is assumed to be of the form X(t) = [q;dq]
            
            qDim = length(X)/2;
            
            % dXdt is of the form dX(t) = [dq;ddq];
            dXdt = [X(qDim+1:2*qDim,1);     % Velocity states
                    U];                     % Acceleration states
        end
        % STATE UPDATE - SINGLE INTEGRATOR
        function [dXdt] = dynamics_singleIntegrator(X,U)
            % The relationship between the input and state vector is by
            % definition of a single integrator: dX/dt = U
            
            % Sanity check
            assert(numel(X) == numel(U),"Expecting the number of inputs to correspond to the state vector.");
            
            % X(t) here is assumed to be of the form X(t) = [q]
            
            % The change in state is literally the input
            dXdt = U;            
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
    end
    
    %% ////////////////// SIMULATION & CORE INTERFACES ////////////////////
    methods
        % INITIAL 3D STATE - [x y z phi theta psi dx dy dz dphi dtheta dpsi]
        function [this] = initialise_3DVelocities(this,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            
            % BUILD STATE VECTOR
            this.localState = zeros(12,1);
            this.localState(4:5) = localXYZrotations(1:2);
            %this.localState(4:6) = localXYZrotations;       % Represent yaw relative to initial position    
            this.localState(7:9) = localXYZVelocity;
            % RETAIN THE PRIOR STATE FOR REFERENCE
            this.VIRTUAL.priorState = this.localState;
        end
        % GLOBAL UPDATE - STATE VECTOR DEFINED AS: [x_t;x_dot]'
        function [this] = updateGlobalProperties_3DVelocities(this,dt,eulerState)
            % This function preforms the state update for a state vector
            % defined as [x y z phi theta psi dx dy dz dphi dtheta dpsi].
            
            positionIndices = 1:3;
            angleIndices = 4:6;
            
            % USE THE 'RATE' STATES DIRECTLY
            localLinearRates  = eulerState(positionIndices+6);
            localAngularRates = eulerState(angleIndices+6);
            
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                                          % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = this.GetVIRTUALparameter('quaternion');   
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
            globalPosition_k_plus = this.VIRTUAL.globalPosition + dt*this.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            this = this.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                    globalVelocity_k_plus,... % Global velocity at k plus
                                                    quaternion_k_plus,...     % Quaternion at k plus
                                                    eulerState);              % The new state for reference
        end
        % GLOBAL UPDATE - EULER 6DOF(3DOF) [x y z phi theta psi],[x y psi]
        function [this] = updateGlobalProperties(this,dt,eulerState)
            % This function computes the position, velocity and 
            % quaternion in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % This function is intended for state vectors of the form:
            % - [x y z phi theta psi]
            % - [x y phi]
                                   
            % NOTATION INDICIES
            if this.Is3D
                positionIndices = 1:3;
                eulerIndices = 4:6;                     % The rotations about X,Y,Z
            else
                positionIndices = 1:2;
                eulerIndices = 3;                       % The rotations about Z
            end
            
            % Populate the prior state if necessary
            if isempty(this.VIRTUAL.priorState)
                this = this.SetVIRTUALparameter('priorState',this.localState);
            end
            % Get the previous state
            priorState = this.GetVIRTUALparameter('priorState');
            
            % Equivalent rates
            localLinearRates  = (eulerState(positionIndices) - priorState(positionIndices))/dt;
            localAngularRates = (eulerState(eulerIndices) - priorState(eulerIndices))/dt;  
            
            % EULER ROTATIONS -> GLOBAL ROTATIONS
            if numel(localAngularRates) ~= 3
                localAngularRates = [0;0;1]*localAngularRates; 
                localLinearRates  = [localLinearRates;0];
            end
                        
            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                   % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = this.GetVIRTUALparameter('quaternion');   
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
            globalPosition_k_plus = this.VIRTUAL.globalPosition + dt*this.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            this = this.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                    globalVelocity_k_plus,... % Global velocity at k plus
                                                    quaternion_k_plus,...     % Quaternion at k plus
                                                    eulerState);              % The new state for reference
        end
        % GLOBAL UPDATE - EULER 6DOF(3DOF) (NO ROTATIONS)
        function [this] = updateGlobalProperties_fixedFrame(this,dt,eulerState)
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
            
            % Equivalent velocity
            velocity_k_plus = (eulerState(positionIndices) - this.VIRTUAL.priorState(positionIndices))/dt;
            
            % ROTATION RATES ABOUT THE GLOBAL AXES
            if numel(eulerState) ~= 6
                velocity_k_plus = [velocity_k_plus;0];
            end
            % Get the current rotation
            q_k = this.GetVIRTUALparameter('quaternion');
            % NEW ROTATION MATRIX (G>B)            
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(q_k);
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = this.VIRTUAL.globalPosition + dt*this.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end
        % GLOBAL UPDATE - DIRECTLY (global default)
        function [this] = updateGlobalProperties_direct(this,p,v,q,X)
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
            assert(numel(p) == 3 && size(p,2) == 1,'Global position must be a 3D column vector [3x1].');
            assert(numel(v) == 3 && size(v,2) == 1,'Global velocity must be a 3D column vector [3x1].');
            assert(numel(q) == 4 && size(q,2) == 1,'Global pose must be a 4D quaternion vector [4x1].');
            assert(numel(X) == numel(this.localState) && size(this.localState,2) == 1,'The length of the objects state update must match the its local state.');       
                                    
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            % Assign the global parameters
            this.VIRTUAL.globalPosition = p;                                % Reassign the global position
            this.VIRTUAL.globalVelocity = v;                                % Reassign the global velocity
            this.VIRTUAL.quaternion = q;                                    % Reassign the quaternion
            this.VIRTUAL.R = OMAS_geometry.quaternionToRotationMatrix(q);   % K_plus rotation
            this.VIRTUAL.priorState = this.localState;                       % Record the previous state
            this.localState = X;                                            % Update the current state
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
                entityRadius = entity.GetVIRTUALparameter('radius');
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
            ind = strfind(repoString,'objects');
            repoString = repoString(1:(ind-1));
            
            % Add the path to the repo's environment 
            addpath([repoString,'environment']);
            
            % Check all the OMAS_dependencies
            if ~OMAS_system.GetFileDependancies()
                error('Unable to get OpenMAS dependancies.');
            end
        end
    end
end