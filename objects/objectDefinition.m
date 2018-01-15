%% GENERIC OBJECT CLASS (objectDefinition.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic object and import its variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite 20/05/2016

classdef objectDefinition
    %% BASIC OBJECT PROPERTIES
    properties
        % SIMULATION IDENTIFIERS
        objectID;
        name;
        mass = 0;
        % QUATERNION UPDATED INDEPENDANTLY OF THE STATE
        globalPosition = zeros(3,1);
        globalVelocity = zeros(3,1);
        quaternion = [1;0;0;0];      % Quaternion - describes the local axis to global rotation
        % KINEMATIC PROPERTIES
        localState  = zeros(12,1);   % Other properties derived from the state.
        % VIRTUAL PROPERTIES
        VIRTUAL;                     % Virtual data container
    end
    %%  CLASS METHODS
    methods (Access = public)
        %% CONSTRUCTOR
        function obj = objectDefinition(varargin)
            % This function is designed to automatically apply a generic
            % naming convention to all objects generated prior to a
            % simulation.
            % INPUTS:
            % varargin - A cell array of property/value pairs 
            % OUTPUTS:
            % obj     - The generated object
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % GENERATE OBJECT-ID (If NOT ALLOCATED BY A SUPERCLASS)
            spawnNewObject = 0;
            if isempty(obj.objectID)
                spawnNewObject = 1;
                [obj.objectID] = objectDefinition.getObjectID();           % Assign the number as an ID tag
            end
            % GENERATE OBJECT NAME
            if isempty(obj.name) || strcmp(obj.name,'')
                % No input name string specified, use greek naming scheme
                [obj.name] = obj.getObjectName(obj.objectID);
            end
                        
            % DEFINE DEFAULT VIRTUAL PROPERTIES
            obj.VIRTUAL.type = OMAS_objectType.misc;                       % Indicate whether the agent is 'active' or 'passive'
            obj.VIRTUAL.size = 1;                                          % Virtual size
            obj.VIRTUAL.colour = rand(1,3);                                % Virtual colour (for plotting)
            obj.VIRTUAL.symbol = 'square';
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
            
            % DISPLAY THE CURRENT AGENT PROPERTIES
            if spawnNewObject
                fprintf('objectcount: %d\tname: %s\ttype: %s\n',obj.objectID,obj.name,class(obj));
            end
        end 
        
        %% INITIALISATION MECHANISMS
        % SET OBJECT LOCAL STATES     [x;y;z;u;v;w;psi;the;phi;p;q;r]
        function [obj] = set_localState(obj,localState)
            % This function aims to intialise an agent entity with a 
            % defined 12 dimensional state. The state is expected in the
            % form [x;y;z;v;u;w;psi;the;phi;p;q;r]. But can also be
            % provided as [x;y;z;v;u;w] for the non-rotational case.
            % INPUT:
            % state - The initial state vector
            %
            % OUTPUT:
            % obj   - The initised agent object with associated properties
            
            % INPUT HANDLING
            if nargin == 0 || ~exist('state','var')
                localState = zeros(12,1);                       % No input given 
            elseif 1 == obj.validateVector(localState,12,1)
                error('State allocation invalid (%s), 12 element state vector required.',obj.name);                   
            end
            % ALLOCATE STATE VARIABLES
            obj.localState = localState;                               % Assign state vector
            obj.quaternion = obj.fromRotations(obj.localState(7:9,1)); % Calculate quaternion from angles
        end    
        % INITIALISE WITH POSITION  [x;y;z]
        function [obj] = set_quaternion(obj,quaternion)
            % INPUT HANDLING
            if nargin == 0 || length(quaternion) ~= 4          
                error('The quaternion is not valid, (%s)',num2str(length(quaternion)));
            end
            % ASSIGN ATTITUDE PROPERTIES
            obj.quaternion = quaternion;                            % Assign quaternion
            obj.localState(7:9,1) = obj.toRotations(quaternion);    % Calculate angles from quaternion
        end
        
        %% UPDATE MECHANISMS
        % THE GENERIC OBJECT CYCLE
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
            linearAcceleration = [0;0;0];
            angularAcceleration = [0;0;0];
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_accelerations(dt,linearAcceleration,angularAcceleration);
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
        
        % BASIC DYNAMICS/KINEMATICS
        % The agent object has both a local NED and global ENU state to
        % allow control and avoidance simulation to be conducted 
        % simulaneously within the common frame of references. The local
        % state vector is of the form: 
        % Xi_frd = [x;y;z;u;v;w;psi;the;phi;p;q;r]_i
        % where as the global is of the form:
        % Xi_enu = [x;y;z;dx;dy;dz;q1;q2;q3;q4]_i.
        
        % The functions below provide the utilities to update both during
        % each time cycle.
        
        % UPDATE AGENT STATE VECTOR FROM LOCAL ACCELERATIONS
        function newState = stateDynamics_accelerations(obj,dt,linearAcceleration,angularAcceleration)
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
                       
            % CHECK DIMENSION OF THE STATE UPDATE
            inputDim = numel(linearAcceleration);
            % DEFINE THE VELOCITY UPDATE
            accelerationVector = zeros(6,1);
            if inputDim < 3                                                % The 2D case
                accelerationVector(6,1) = angularAcceleration;             % is [1x1]
            else                                                           % The 3D case
                accelerationVector(4:(3+inputDim),1) = angularAcceleration;% is [3x1]
            end
            accelerationVector(1:inputDim,1) = linearAcceleration;         % Linear velocity
            % COMPUTE CHANGE IN STATES
            dtMatrix = dt*eye(12);
            % STATE VECTOR = CONCATINATED ZEROTH AND FIRST ORDER TERMS
            dX = [obj.localState(7:12,1);accelerationVector];
            newState = obj.localState + dtMatrix*dX;                   % Access the current state and append difference
        end
        % UPDATE AGENT STATE VECTOR FROM LOCAL VELOCITIES
        function newState = stateDynamics_velocities(obj,dt,linearVelocity,angularVelocity)
            % This function provides a basic model for the agent
            % kinematics, controlled by linear and angular velocity inputs
            % INPUTS:
            % dt              - The simualation timestep
            % linearVelocity  - The agent's current linear velocity [2x1] or [3x1]
            % angularVelocity - The agent's current angular velocity [1x1] or [3x1]
            % OUTPUTS:
            % newState        - The updated local (FRD) state vector [12x1]
            
            % CHECK DIMENSION OF THE STATE UPDATE
            inputDim = numel(linearVelocity);
            % DEFINE THE VELOCITY UPDATE
            velocityVector = zeros(6,1);
            if inputDim < 3                                                % The 2D case
                velocityVector(6,1) = angularVelocity;                     % is [1x1]
            else                                                           % The 3D case
                velocityVector(4:(3+inputDim),1) = angularVelocity;        % is [3x1]
            end
            velocityVector(1:inputDim,1) = linearVelocity;                 % Linear velocity
            % GET THE PRIOR STATE
            priorState = obj.localState;
            priorState(7:12) = velocityVector;
            % INTEGRATATION MATRIX
            dtMatrix = diag(ones(size(obj.localState)));
            dtMatrix(1:6,7:12) = dt*eye(6); 
            % INTEGRATE THE PREVIOUS STATE AND THE CHANGE (X = X + dX)
            newState = dtMatrix*priorState;
        end
        % REDEFINE THE OBJECTS CURRENT GLOBAL PROPERTIES
        function [obj] = updateGlobalProperties(obj,dt,localFRDState)
            % This function calculates an objects equivalent global
            % representation (position, velocity and quaternion) based on 
            % its new FORWARD-RIGHT-DOWN (LOCAL NED) state vector: 
            % [x; y; z; phi; theta; psi; u; v; w; phi_dot; theta_dot; psi_dot]
            
            % INPUT VALIDATION
            if ~isreal(localFRDState) || ~isnumeric(localFRDState)
                error('The resulting state must be real, or is invalid');
            end
            
            % REASSIGN NEW LOCAL CONTROL (NED) STATE
            obj.localState = localFRDState;
            % DEFINE LOCAL STATE INDICES
            positionIndices = 1:3;
            rotationIndices = 4:6;
            velocityIndices = 7:9;
            omegaIndices = 10:12;
            % WE MUST DEFINE THE LOCAL FORWARDS-LEFT-UP REPRESENTATION 
            localFLUState = localFRDState;
            % LOCAL POSITIONS
            localFLUState(positionIndices(2)) = -localFLUState(positionIndices(2)); % Negate Y NED position
            localFLUState(positionIndices(3)) = -localFLUState(positionIndices(3)); % Negate Z NED position
            % LOCAL ROTATION
            localFLUState(rotationIndices(1)) = localFRDState(rotationIndices(1)) - pi; % ROTATE NED INTO FLU (180 IN ROLL)
            % LOCAL VELOCITIES
            localFLUState(velocityIndices(2)) = -localFLUState(velocityIndices(2)); % Negate Y NED velocity
            localFLUState(velocityIndices(3)) = -localFLUState(velocityIndices(3)); % Negate Z NED velocity 
            % LOCAL RATES
            localFLUState(omegaIndices(2)) = localFLUState(omegaIndices(2));  % Negate rotation about y axis
            localFLUState(omegaIndices(3)) = -localFLUState(omegaIndices(3)); % Negate rotation about z axis            
            % CONVERT THE FLUE STATES TO UP
            XYZrates(1) = localFLUState(omegaIndices(3));    % Map rotations around rotated body Z axis to global Z
            XYZrates(2) = localFLUState(omegaIndices(2));    % Map rotations around rotated body Y axis to global Y
            XYZrates(3) = localFLUState(omegaIndices(1));    % Map rotations around rotated body X axis to global X
            % UPDATE THE GLOBAL QUATERNION
            [obj.quaternion] = OMAS_axisTools.updateQuaternion(obj.quaternion,XYZrates,dt);
            % GET UPDATED INVERSE ROTATION MATRIX
            [R_BG,~] = OMAS_axisTools.quaternionToRotationMatrix(obj.quaternion);
            % UPDATE GLOBAL POSE
            velocityUpdate = R_BG*localFLUState(velocityIndices,1);
            obj.globalPosition = obj.globalPosition + velocityUpdate*dt;
            obj.globalVelocity = velocityUpdate;
        end
        
        % OBJECT CONFIGURATION PARSER
        function [obj] = configurationParser(obj,objectParameters)
            % This function is designed to parse a generic list of input
            % pararmeters and compare them to properties of the current
            % object.
                        
            % If a valid numeric vector is provided; substitute
            if exist('objectParameters','var') && iscell(objectParameters) % If the first input is a structure
                virtualFlag = isprop(obj,'VIRTUAL');                       % Object has virtual properties
                
                for parameterIndex = 1:numel(objectParameters)
                    % FOR EACH USER-DEFINED PARAMETER
                    givenProperty = objectParameters{parameterIndex};
                    if ~ischar(givenProperty)
                        continue
                    else
                        % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                        if isprop(obj,givenProperty)
                            obj.(givenProperty) = objectParameters{parameterIndex + 1};   % Make a substitution
                        end
                        % IF THERE EXISTS A VIRTUAL PROPERTY BY THAT NAME
                        if virtualFlag && isfield(obj.VIRTUAL,givenProperty) 
                            obj.VIRTUAL.(givenProperty) = objectParameters{parameterIndex + 1};
                        end
                    end
                end
            end
        end
    end
    % STATIC METHODS & TOOLS
    methods (Static)
        %% GENERAL USER TOOLS
        % BOUND VALUE
        function [boundedVector] = boundValue(inputVector,lower,upper)
            % This function simply bounds a value between two given values
            
            stateNumber = size(inputVector,1);
            
            if stateNumber ~= size(lower,1) || stateNumber ~= size(upper,1)
                error('Input vector must have bounds for each dimension')
            end
            
            % ELSE ASSUME THERE IS UPPER AND LOWER BOUND FOR EACH VECTOR
            % DIMENSION
            for ind = 1:stateNumber
                if inputVector(ind) < lower(ind)
                    inputVector(ind) = lower(ind);
                end
                if inputVector(ind) > upper(ind)
                    inputVector(ind) = upper(ind);
                end
            end
            boundedVector = inputVector;
        end
    
        %% GENERAL ROTATION & MATH TOOLS
        % GET QUATERNION UPDATE FROM AXIS RATES
        function [q] = updateQuaternion(q0,axisRates,dt)
            % This function updates a given quaternion using the body axis
            % rates.
            % INPUTS:
            % q0        - The initial quaternion
            % axisRates - The body axis rates
            % dt        - The time displacement
            
            % Get the normalising diagonal term
            correction = 1 - sqrt(sum(q0.^2));
            % Calculate the quaternion difference
            q_dot = 0.5*[correction,  axisRates(3),  -axisRates(2),  axisRates(1);
                      -axisRates(3),    correction,   axisRates(1),  axisRates(2);
                       axisRates(2), -axisRates(1),     correction,  axisRates(3);
                      -axisRates(1), -axisRates(2),   -axisRates(3),   correction]*q0;
            % Integrate the quaternion difference
            q = q_dot*dt + q0;
            % Re-normalise
            q = q/sqrt(sum(q.^2));
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R_AB,R_BA] = rotationMatrix(q)
            % This function defines the equivalent rotation matrix from the
            % provided quaternion. (Associated block "Quaternion to DCM")
            % INPUT:
            % q    - The quaternion rotation
            % OUTPUT:
            % R    - The rotation matrix through the same quaternion
            % Rinv - The rotation matrix describing the reverse rotation
            
            R_AB = zeros(3,3);           
            % Normalise the quaternion
            q = OMAS_axisTools.qUnit(q);
            % Assemble the quaternion elements
            R_AB(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
            R_AB(1,2) = 2*(q(1)*q(4) + q(2)*q(3));
            R_AB(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
            R_AB(2,1) = 2*(q(3)*q(2) - q(1)*q(4));
            R_AB(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
            R_AB(2,3) = 2*(q(1)*q(2) + q(3)*q(4));
            R_AB(3,1) = 2*(q(1)*q(3) + q(2)*q(4));
            R_AB(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
            R_AB(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
            % Transpose the matrix for the global-body transformations
            R_BA = transpose(R_AB);
        end
        % CONVERT ROTATION ANGLES INTO A QUATERNION    
        function [q] = fromRotations(rotationAngles)
            % THis function converts euler angle rotations into a
            % quaternion via its components. Associated block:
            % "Rotation Angles to Quaternions".
            
            rotationAngles = 0.5*rotationAngles;
            % Build vector of trigonometric arguements
            trigArgs = [sin(rotationAngles);cos(rotationAngles)];
            % Assemble Quaternion components
            q = zeros(4,1);
            q(1) = trigArgs(4)*trigArgs(5)*trigArgs(6) + trigArgs(1)*trigArgs(2)*trigArgs(3);
            q(2) = trigArgs(4)*trigArgs(5)*trigArgs(3) - trigArgs(1)*trigArgs(2)*trigArgs(6);
            q(3) = trigArgs(4)*trigArgs(2)*trigArgs(6) + trigArgs(1)*trigArgs(5)*trigArgs(3);
            q(4) = trigArgs(1)*trigArgs(5)*trigArgs(6) - trigArgs(4)*trigArgs(2)*trigArgs(3);
        end
        % CONVERT QUATERNION INTO ROTATION ANGLES
        function [angles] = toRotations(q)
            % Gets the rotation angles from an equivalent quaternion
            % attitude. Associated block:
            % "Quaternions to Rotation Angles"
            if length(q) ~= 4
                error('Input quaternion is of length %s',num2str(length(q)));
            end
            
            % Normalise the quaternion
            q = quaternionTools.unit(q);
            angles = zeros(3,1);
            % Define the euler angles
            angles(1) = atan2((2*(q(2)*q(3) + q(1)*q(4))),(q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2));  % Phi   (roll)
            angles(2) = asin(-2*(q(2)*q(4) - q(1)*q(3)));                                       % Theta (pitch)
            angles(3) = atan2((2*(q(3)*q(4) + q(1)*q(2))),(q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2));  % Psi   (yaw)
        end
        % ROTATE VECTOR THROUGH AN ANGLE, AROUND A GIVEN AXIS
        function [newVector] = rotateVectorAboutAxis(oldVector,axisVector,theta)
            % This function is designed to calculate a vector
            % following a rotation around a given axis vector, through a
            % given angle.
            % INPUTS:
            % oldVector  - The initial vector
            % axisVector - The axis of rotation
            % theta      - The angle of rotation
            % OUTPUTS:
            % newVector  - The rotated 3D vector
            
            % NORMALISE THE AXIS VECTOR
            axisVector = axisVector/norm(axisVector);  % Normalize rotation axis
            % GET THE CROSS PRODUCT PROECTION BETWEEN THE AXIS AND VECTOR
            crossVector = cross(axisVector,oldVector);
            
            % DETERMINE THE MAPPING OF EACH VECTOR COMPONENTS
            newVector = cos(theta)*oldVector ...
                + (crossVector)*sin(theta)  ...
                + axisVector*(dot(axisVector,oldVector))*(1 - cos(theta));
        end
        % GET THE TENSOR PRODUCT OF TWO n-DIMENSIONAL VECTORS
        function Z = tensorProduct(X,Y)
            %   Z = TensorProduct(X,Y) returns the REAL Kronecker tensor product of X and Y.
            %   The result is a multidimensional array formed by taking all possible products
            %   between the elements of X and those of Y.
            %
            %   If X is m-by-n and Y is p-by-q-by-r, then kron(X,Y)
            %   is m-by-p-by-n-by-q-by-r.
            %
            %   X and Y are multidimensional matrices
            %   of any size and number of dimensions
            %
            %   E.g. if X is of dimensions (4, 5, 3) and Y of dimension (3, 1, 7, 4)
            %   TensorProduct(X, Y) returns a multidimensional matrix Z of dimensions:
            %   (4, 5, 3, 3, 7, 4)
            %
            %   $ Date: 2001/11/09 10:20:00 GMT $
            %
            %       Steeve AMBROISE --> sambroise@gmail.com
            %
            sX=size(X);sY=size(Y);
            ndim1=ndims(X);ndim2=ndims(Y);
            indperm=[ndim2+1:ndim1+ndim2,1:ndim2];
            % to remove all singleton dimensions
            Z=squeeze(repmat(X,[ones(1,ndims(X)),sY]).*...
                      permute(repmat(Y,[ones(1,ndims(Y)),sX]),indperm));
        end
        % CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
        function [outputMatrix] = skew(inputVector)
            % This function generates a skew-symmetric for the computation
            % of the vector cross-product.
            % INPUTS:
            % inputVector - The original 3D vector
            % OUTPUT:
            % outputMatrix - The equivalent skew-symmetric matrix
            
            if length(inputVector) ~= 3
                warning('The input vector must be three dimensional.');
                return
            end
            % Apply element mapping
            outputMatrix = zeros(3,3);
            outputMatrix(1,2) = -inputVector(3);
            outputMatrix(1,3) =  inputVector(2);
            outputMatrix(2,1) =  inputVector(3);
            outputMatrix(2,3) = -inputVector(1);
            outputMatrix(3,1) = -inputVector(2);
            outputMatrix(3,2) =  inputVector(1); % Arrange the components
        end
        % UNIT VECTOR
        function [unitVector] = unit(inputVector)
            % This function returns the unit vector of a vector
            unitVector = inputVector/sqrt(sum(inputVector.^2));
        end
    end
    % SIMULATION UTILITIES 
    methods (Static, Access = private)
        % SIMULATION TOOLS
        % APPLY NAMING CONVENTION
        function [namestr] = getObjectName(objectID)
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
    % PRIVATE METHODS (CLASS TOOLS)
    methods (Access = private)        
        % GENERATE THE OBJECT FIGURE HANDLE
        function [objectSurface] = getObjectFigureHandle(obj,state)  
            %% Determine if a figure is already open:
            % BUILD REFERENCE OBJECTS
            [X,Y,Z] = sphere(15);
            X = X.*obj.VIRTUAL.size + state(1,1);
            Y = Y.*obj.VIRTUAL.size + state(2,1);
            Z = Z.*obj.VIRTUAL.size + state(3,1);
            % Append xyz grid data
            objectSurface = struct('X',X,'Y',Y,'Z',Z);
        end
    end
end

