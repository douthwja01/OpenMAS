%% SCENARIO BUILDER CLASS (scenarioBuilder.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The scenario builder is a tool designed to provide a set of tools to aid
% in the design of object scenarios in the global scene. 

classdef scenarioBuilder
    % BASIC SCENARIO DEFINITION
    % This class handles the methods of generating different scenarios of
    % a defined number of objects/agents.
    
    % SCENARIO PROPERTIES
    properties
        % META DATA
        objects;
        name;
        % GLOBAL DATA
        % ALL GLOBAL DATA IS GIVEN IN THE EAST-NORTH-UP (ENU) FRAME OF
        % REFERENCE
        globalTriad = [1,0,0;0,1,0;0,0,1];      % XYZ (Earth) frame of reference
        position;                               % Global positions
        velocity;                               % Global velocities
        quaternion;                             % The quaternion moving from the global XYZ
    end
    % PUBLIC METHODS
    methods  
        % CONSTRUCTION METHOD 
        function obj = scenarioBuilder(objectNumber)
            % Default scenario with one agent
            if nargin == 0 || 0 == isnumeric(objectNumber)
                objectNumber = 1; % Provide a
            end
            % DECLARE ANY INITIAL PARAMETERS
            obj.name = 'Default Scenario.';
            obj.objects = objectNumber;
            
            % POPULATE PLACEHOLDER VALUES
            obj.position   = zeros(3,objectNumber);
            obj.velocity   = zeros(3,objectNumber);
            obj.quaternion = repmat([1;0;0;0],[1 objectNumber]);
        end
               
        %% SCENARIO GENERATORS
        % PLANAR CONCENTRIC DISK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = planarDisk(obj,varargin)
            % This function generates a configuration representing a given 
            % number of objects distributed about a disk defined by a 
            % central axis vector.
            % INPUTS:
            % pointA      - The first axis reference point.
            % pointB      - The second axis reference point; the ring center.
            % radius      - The scalar radius of the ring
            % velocity    - The velocity multiplier for the object set.
            % OUTPUTS:
            % definition  - The scenario definition
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct('pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'velocity',0,...
                                   'scale',1);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % ////////////// DESIGN THE DISK DISTRIBUTION /////////////////
            objectNumber = obj.objects;
            % If there is only one object, simply place it at the center
            if objectNumber == 1
                scenarioConfig.position(:,1) = config.pointB;
                scenarioConfig.velocity(:,1) = [1;0;0]*config.velocity;
                scenarioConfig.quaternion(:,1) = [1;0;0;0];
                return
            end
            % The configuration container for a disk of objects.            
            scenarioConfig = struct('objects',objectNumber,...
                                'name','Distributed planar disk.',...
                                'position',[],...
                                'velocity',[],...
                                'quaternion',[]);
            remainingObjects = objectNumber;
            layer = 1; 
            while layer < objectNumber && remainingObjects > 0
                % NEW LAYER RADIUS
                layerRadius = layer*config.scale;                          % New circle radius
                criticalLayerAngle = 2*asin(1/(2*layer));                  % New maximum nodal angle  
                layerMax = floor(2*pi/criticalLayerAngle);                 % Maximum number on the new layer
                % DISTRIBUTE OBJECTS OVER THE LAYERS
                if layerMax > remainingObjects
                    objectsInLayer = remainingObjects;
                else
                    objectsInLayer = layerMax;
                end
                % BUILD THE NEXT LAYER
                [~,layerConfig] = obj.planarRing('objects',objectsInLayer,...
                                                 'pointA',config.pointA,...
                                                 'pointB',config.pointB,...
                                                 'zeroAngle',config.zeroAngle,...
                                                 'radius',layerRadius);
                % EXTRACT THE LAYER GLOBAL PROPERTIES
                first = 1 + objectNumber - remainingObjects;
                last  = (first - 1) + objectsInLayer;
                indexVector = first:last;
                scenarioConfig.position(:,indexVector) = layerConfig.position(:,1:objectsInLayer);
                scenarioConfig.velocity(:,indexVector) = layerConfig.velocity(:,1:objectsInLayer);
                scenarioConfig.quaternion(:,indexVector) = layerConfig.quaternion(:,1:objectsInLayer);
                % INCREMENT LAYER
                remainingObjects = remainingObjects - objectsInLayer;
                layer = layer + 1;
            end

        end
        % PLANAR RING OF STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = planarRing(obj,varargin)
            % This function generates a configuration representing a given 
            % number of objects destributed around a ring defined by a 
            % central axis vector. 
            % INPUTS:
            % pointA     - The first axis reference point.
            % pointB     - The second axis reference point; the ring center.
            % mod_radius - The scalar radius of the ring
            % velocity   - The velocity multiplier for the object set.
            % OUTPUTS:
            % definition - The scenario definition
                        
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',obj.objects,...
                                   'pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'radius',10,...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'offsetAngle',pi/5,...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);

            % DEFINE AN ANGLE FOR EQUAL DISTRIBUTION ABOUT THE RING
            nodalAngle = (2*pi)/config.objects;          
            % HAND NODAL ANGLE TO ringAngle FUNCTION
            [ obj, scenarioConfig ] = obj.planarAngle('objects',config.objects,...
                                                      'pointA',config.pointA,...
                                                      'pointB',config.pointB,...
                                                      'radius',config.radius,...
                                                      'zeroAngle',config.zeroAngle,...
                                                      'offsetAngle',nodalAngle,...
                                                      'velocity',config.velocity);
            % UPDATE SCENARIO LABEL
            obj.name = 'Planar ring.';
        end
        % CONSTANT OFFSET ANGLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = planarAngle(obj,varargin)
            % GENERATE A SCENARIO DEFINITION 
            % INPUTS:
            % pointA        - The first axis reference point.
            % pointB        - The second axis reference point; the ring center.
            % mod_radius    - The scalar radius of the ring
            % velMultiplier - The velocity multiplier for the object set.
            % offsetAngle   - The constant radial angle 
            % OUTPUTS:
            % definition    - The scenario deifnition
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',obj.objects,...
                                   'pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'radius',10,...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'offsetAngle',pi/5,...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            % UPDATE SCENARIO LABEL
            obj.name = sprintf('Planar offset angle of %srad.',num2str(config.offsetAngle));
            
            % GENERATE THE NODE(OBJECT) CARTAESIAN POSITIONS      
            % GET THE AXIS VECTOR PROPERTIES
            localXAxis = config.pointB - config.pointA;
            unit_axisVector = OMAS_geometry.unit(localXAxis);
            % GET THE RADIAL VECTOR
            perpVector = cross(unit_axisVector,[1;0;0]);
            if sum(perpVector) == 0
                perpVector = cross(unit_axisVector,[0;1;0]);
            end
            perpVector = OMAS_geometry.unit(perpVector); % Re-normalise the perpendicular vector
            % SCALE THE UNIT RADIAL VECTOR
            perpVector = config.radius*perpVector;               
            
            % GET THE NODE POINT SET
            nodalAngle = config.zeroAngle + pi/2; % pi/2 aligns the first agent with the x-axis
            for node = 1:config.objects
                % ROTATE THE RADIAL VECTOR TO DEFINE NODAL POSITIONS
                radialVector = OMAS_geometry.rodriguesRotation(perpVector,unit_axisVector,nodalAngle);
                % NODAL POSITIONS
                globalPosition = config.pointB + radialVector;
                % NODAL VELOCITIES
                unitVelocity = - OMAS_geometry.unit(radialVector);        % -ve sign to make concentric velocities +ve             
                % NODAL ORIENTATIONS
                % We desire to orientate each object towards the center of
                % the ring. This is done by aligning the local body XY plane
                % With the plane of the ring. The local X axis is then
                % aligned with the cocentric global velocity. 
                                                
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                localTriad = zeros(3);
                localTriad(:,1) = OMAS_geometry.unit(unitVelocity);             % GET THE UNIT DIRECTION VECTOR
                localTriad(:,3) = OMAS_geometry.unit(unit_axisVector);          % UNIT RING AXIS
                localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));
                
                % STORE THE OBJECTS IN SCENARIO DEFINITION 
                obj.position(:,node) = globalPosition;
                obj.velocity(:,node) = config.velocity*unitVelocity; 
                obj.quaternion(:,node) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                
%                 % ROBOTICS TOOLBOX: VALIDATION
%                 Rxyz = quat2rotm(transpose(obj.quaternion(:,node)));
%                 eulA = rotm2eul(Rxyz,'XYZ');
%                 % OpenMAS TOOLBOX: TO VALIDATE
%                 Rxb = OMAS_geometry.quaternionToRotationMatrix(obj.quaternion(:,node));
%                 eulB = OMAS_geometry.rotationMatrixToEulers(Rxb);

                % Object data is concatinated vertically
                % Increment the nodal position by the offset angle
                nodalAngle = nodalAngle - config.offsetAngle;
            end

            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',config.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...     % The global position
                                    'velocity',obj.velocity,...     % The velocity of the object in the global frame
                                    'quaternion',obj.quaternion);   % The attitude of the body in the global axes
            % UPDATE SCENARIO
            obj.name = 'Planar Angle.';
        end
        % SPHERICAL EQUAL DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = regularSphere(obj,varargin)
            % This function is designed to generate points equally
            % distributed about the circumference of a sphere.
                       
            % DEFAULT CONFIGURATION
            defaultConfig = struct('radius',10,...
                                   'center',[0;0;0],...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
                        
            % ATTEMPT TO ADD THE ICOSAHEDRALS LIBRARY
            try
                libraryPath = strcat(cd,'\scenarios\icosahedrals');
                addpath(libraryPath);
            catch
                error('Unable to load the icosahedral generator function, has it been moved?');
            end                 
            
            % GET SCALED ICOSAHEDRAL SPHERE
            % This function will generate a vector of points equally
            % distributed about the circumference of a sphere.
            [pointCloud,~] = spheretri(obj.objects);
            XYZ = (pointCloud.*config.radius)';              
            XYZ = XYZ + config.center;                                      % Offset and scale the sphere
            
            % SORT THE POINTS
            XYZ = unique(XYZ','rows');                                     % Remove repeat entries
            XYZ = sortrows(XYZ,1);
            XYZ = XYZ';
%             planarIndices = any(XYZ == 0);
%             nodalPositions = XYZ(:,planarIndices);
            nodalPositions = XYZ();
            nodalPositions = nodalPositions(:,1:obj.objects);
            % ALLOCATE RANDOM SET TO OBJECT GLOBAL POSITIONS
            for index = 1:obj.objects               
                % DETERMINE THE NODAL POSITIONS
                nodalPosition = config.center + nodalPositions(:,index);
                
                errorFlag = 1;
                while errorFlag == 1   
                    % DEFINE GLOBAL VELOCITIES
                    unit_radii = OMAS_geometry.unit(config.center - nodalPosition);
                    obj.velocity(:,index) = config.velocity*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = OMAS_geometry.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.position(:,index) = nodalPosition;
                        % ROTATION
                        obj.quaternion(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                        errorFlag = 0;
                    catch
                        % IF ROTATION IS INVALID, RESAMPLE THE POSITION
                        [nodalPosition,~] = datasample(XYZ,1,2,'Replace',false);
                    end
                end
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);            
            % UPDATE SCENARIO LABEL
            obj.name = 'Equally spaced spherical distribution.';   
        end
        % COCENTRIC HELICAL STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = helix(obj,varargin)
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct('file','scenario.mat',...
                                   'objects',obj.objects,...
                                   'pointA',[0;0;0],...
                                   'pointB',[0;0;10],...
                                   'radius',10,...
                                   'offsetAngle',pi/5,...
                                   'waypointRadius',2,...
                                   'velocities',0,...
                                   'plot',0);
            % PARSE CONFIGURATIONS
            [config] = scenarioBuilder.configurationParser(defaultConfig,varargin);
            
            % DEFINE THE HELIX AXIS
            helixAxis = config.pointB - config.pointA;
            norm_helixAxis = norm(helixAxis);
            unit_helixAxis = helixAxis/norm_helixAxis;
            
            % THE UNIT SPACING
            objectNumber = obj.objects;
            axialSeparation = norm_helixAxis/objectNumber;
            
            % REFERENCE POSITION FOR FIRST NODE
            pointA_i = config.pointA - axialSeparation*unit_helixAxis;
            for i = 1:objectNumber
                % AXIAL PARAMETERS
                pointB_i = pointA_i + axialSeparation*unit_helixAxis;
                radialAngle = config.offsetAngle*(i-1);
                
                % DEFINE CONFIGURATION FOR OBJECT "i"
                [ ~,planarConfig ] = obj.planarAngle('objects',1,...                   
                                                     'pointA',pointA_i,...
                                                     'pointB',pointB_i,...
                                                     'radius',config.radius,...
                                                     'zeroAngle',radialAngle,...                       % The rotation of the first node
                                                     'offsetAngle',0);   
                % DEFINE THE CONFIGURATION OF THE AGENT
                obj.position(:,i) = planarConfig.position(:,1);
                obj.velocity(:,i) = planarConfig.velocity(:,1);
                obj.quaternion(:,i) = planarConfig.quaternion(:,1);
                % SHIFT UP THE AXIS VECTOR
                pointA_i = pointB_i;    
            end
            
            obj.name = 'Cocentric helix.'; 
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion); 
        end
        % COLINEAR DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = line(obj,varargin)
            % DEFAULT CONFIGURATION
            defaultConfig = struct('file','scenario.mat',...
                                   'objects',obj.objects,...
                                   'pointA',[0;0;0],...
                                   'pointB',[0;0;10],...
                                   'heading',[1;0;0],...
                                   'velocities',1,...
                                   'plot',0);
            % PARSE CONFIGURATIONS
            [config] = scenarioBuilder.configurationParser(defaultConfig,varargin);
            
            % LINE AXIS
            axis = config.pointB - config.pointA;
            unit_axis = OMAS_geometry.unit(axis);                          % Get the axis vector between the two points
            
            if sum(iszero(config.heading)) == 3
                headingReference = [1;0;0];
            else
                headingReference = OMAS_geometry.unit(config.heading);
            end
            
            % DISTRIBUTE THE OBJECTS EQUALLY ALONG THE AXIS
            if isnumeric(config.objects)
                nObjects = config.objects;
            else
                nObjects = numel(config.objects);
            end
            % MOVE ALONG THE AXIS AND CREATE OBJECT STATES
            increments = linspace(0,norm(axis),nObjects);
            for node = 1:nObjects
                % OBJECT POSITION
                obj.position(:,node) = config.pointA + increments(node)*unit_axis;
                % OBJECT VELOCITY
                obj.velocity(:,node) = config.velocities;
                % ORIENTATION
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)                  
                localTriad = zeros(3);
                localTriad(:,1) = headingReference;                        % GET THE UNIT DIRECTION VECTOR
                localTriad(:,3) = OMAS_geometry.unit(obj.globalTriad(:,3)); % GLOBAL Z-AXIS
                localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));
                % RESOLVE ROTATION QUATERNION
                obj.quaternion(:,node) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',config.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...     % The global position
                                    'velocity',obj.velocity,...     % The velocity of the object in the global frame
                                    'quaternion',obj.quaternion);   % The attitude of the body in the global axes
        end
        
        % RANDOM SPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = randomSphere(obj,varargin)
            % This fuction is designed to generate random positions 
            % distributed about the circumference of a sphere. 
                        
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',obj.objects,....
                                   'radius',10,...
                                   'center',[0;0;0],...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            % ATTEMPT TO ADD THE ICOSAHEDRALS LIBRARY
            try
                libraryPath = strcat(cd,'\scenarios\icosahedrals');
                addpath(libraryPath);
            catch
                error('Unable to load the icosahedral generator function, has it been moved?');
            end            
            
            % GET SCALED ICOSAHEDRAL SPHERE
            % This function will generate a vector of points equally
            % distributed about the circumference of a sphere.
            [pointCloud,~] = spheretri(obj.objects);
            XYZ = (pointCloud.*config.radius)';              
            XYZ = XYZ + config.center;                                     % Offset and scale the sphere
            XYZ = unique(XYZ','rows');                                     % Remove repeat entries
            XYZ = XYZ';
            
            % RANDOMLY SAMPLE FROM THE POINT SET
            [nodalPositions,~] = datasample(XYZ,obj.objects,2,'Replace',false);   % Take one sample column from X coordinates

            % ALLOCATE RANDOM SET TO OBJECT GLOBAL POSITIONS
            for index = 1:config.objects               
                % DETERMINE THE NODAL POSITIONS
                nodalPosition = config.center + nodalPositions(:,index);
                
                errorFlag = 1;
                while errorFlag == 1   
                    % DEFINE GLOBAL VELOCITIES
                    unit_radii = OMAS_geometry.unit(sphereCenter - nodalPosition);
                    obj.velocity(:,index) = config.velocity*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = OMAS_geometry.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.position(:,index) = nodalPosition;
                        % ROTATION
                        obj.quaternion(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                        errorFlag = 0;
                    catch
                        % IF ROTATION IS INVALID, RESAMPLE THE POSITION
                        [nodalPosition,~] = datasample(XYZ,1,2,'Replace',false);
                    end
                end
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);   
            % UPDATE SCENARIO LABEL
            obj.name = 'Random cocentric sphere';
        end  
        % DEFAULT RANDOM STATE GENERATOR - NORMAL %%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = random(obj)
            % This function is a simple placeholder to facilitate the use
            % of .random for quick access to normally distributed state sets.
            [ obj,scenarioConfig ] = randomNormal(obj);
        end
        % RANDOM (NORMAL DISTRIBUTION) SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = randomNormal(obj)
            % This function generates a random state vector set using a normal distribution. 
            % assuming the order of the states is [x;y;z;u;v;w;phi;theta;psi;p;q;r]
            % INPUT:
            % obj          - The scenario definition object
            % OUTPUT:
            % randomStates - The state vector set [states x objects]
            % obj          - The returned scenario object
            
            % INPUT HANDLING
            if ~isnumeric(obj.objects)
                warning('A number of states number be entered.');
                scenarioConfig = [];
                return
            end
            % DEFINE A RANDOM STATE PARAMETERS
            obj.name = 'Normal-randomised';
            stateNum = 9;
            positionTuple = 1:3;
            velocityTuple = 4:6;
            angleTuple = 7:9;
            stateGains = [10,10,pi];
            normalObjectStates = zeros(stateNum,1);
            
            for index = 1:obj.objects
                % GENERATE RANDOM NUMBER SET
                maxR = max(randn(stateNum,1));
                normalObjectStates(:,index) = randn(stateNum,1)/maxR; % returns an n-by-n matrix of random numbers.
                
                % ALLOCATE STATE VALUES
                normalObjectStates(positionTuple,index) = stateGains(1)*normalObjectStates(positionTuple,index);
                normalObjectStates(velocityTuple,index) = stateGains(2)*normalObjectStates(velocityTuple,index);
                
                % NODAL POSITIONS
                obj.position(:,index) = normalObjectStates(positionTuple,index);
                % NODAL VELOCITIES
                obj.velocity(:,index) = normalObjectStates(velocityTuple,index);
                                                    
                % NODAL ORIENTATIONS
                % We desire to orientate each object towards the center of
                % the ring. This is done by aligning the local body XY plane
                % With the plane of the ring. The local X axis is then
                % aligned with the cocentric global velocity. 
                
                % ALLOCATE ANGULAR STATE VALUES                
                normalObjectStates(angleTuple,index) = stateGains(3)*normalObjectStates(angleTuple,index);
                normalObjectStates(angleTuple(2:end),index) = 0;           % Remove pitch and yaw from random states (for stability)        
                 
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                localTriad = zeros(3);
                localTriad(:,1) = OMAS_geometry.unit(normalObjectStates(velocityTuple,index));      % The body axis X direction
                localTriad(:,3) = OMAS_geometry.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternion(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);
            save('scenario.mat','scenarioConfig');
        end
        % RANDOM (UNIFORM DISTRIBUTION) SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = randomUniform(obj)
            % This function generates a random state vector set using a uniform distribution. 
            % assuming the order of the states is [x;y;z;u;v;w;phi;theta;psi;p;q;r]
            % INPUT:
            % obj          - The scenario definition object
            % OUTPUT:
            % obj          - The returned scenario object
            
            % INPUT HANDLING
            if ~isnumeric(obj.objects)
                warning('A number of states number be entered.');
                return
            end
            % DEFINE A RANDOM STATE PARAMETERS
            obj.name = 'Uniformly-randomised';
            stateNum = 9;
            positionTuple = 1:3;
            velocityTuple = 4:6;
            
            angleTuple = 7:9;
            stateGains = [10,10,pi];
            uniformObjectStates = zeros(stateNum,1);
            
            for index = 1:obj.objects
                % GENERATE RANDOM NUMBER SET
                uniformObjectStates(:,index) = (rand([stateNum,1])-0.5).*2;       % returns an n-by-n matrix of random numbers.
                % ALLOCATE LINEAR STATE VALUES
                uniformObjectStates(positionTuple,index) = stateGains(1)*uniformObjectStates(positionTuple,index);
                uniformObjectStates(velocityTuple,index) = stateGains(2)*uniformObjectStates(velocityTuple,index);
                
                % NODAL POSITIONS
                obj.position(:,index) = uniformObjectStates(positionTuple,index);
                % NODAL VELOCITIES
                obj.velocity(:,index) = uniformObjectStates(velocityTuple,index);
                                
                % NODAL ORIENTATIONS
                % We desire to orientate each object towards the center of
                % the ring. This is done by aligning the local body XY plane
                % With the plane of the ring. The local X axis is then
                % aligned with the cocentric global velocity. 
                
                % ALLOCATE ANGULAR STATE VALUES
                uniformObjectStates(angleTuple,index) = stateGains(3)*uniformObjectStates(angleTuple,index);
                uniformObjectStates(angleTuple(2:end),index) = 0; % Remove pitch and yaw from random states (for stability)
                 
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                localTriad = zeros(3);
                localTriad(:,1) = OMAS_geometry.unit(uniformObjectStates(velocityTuple,index));     % The body axis X direction
                localTriad(:,3) = OMAS_geometry.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = OMAS_geometry.unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternion(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);
            save('scenario.mat','scenarioConfig');
        end
    end

    
    %% GENERIC SCENARIO TOOLS /////////////////////////////////////////
    methods
        % PLOT SCENARIO FROM GLOBAL CONFIGURATION ONLY
        function [figureHandle] = plot(obj)
            % This function is designed to plot the positional and velocity
            % configuration of the current scenario
            % INPUTS:
            % obj         - The scenario object
            % .position   - The position vector set
            % .velocity   - The velocity vector set
            % .quaternion - The quaternion vector set
            
            global plotnum
            
            % DETERMINE PLOT PROPERTIES
            if ~exist('plotnum','var') || isempty(plotnum)
                plotnum = 1;    % Default to first plot
            end
            
            % GENERATE THE FIGURE
            figureHandle = figure(plotnum);
            title('Scenario Definition');
            axis('equal');
            xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
            view([70 25]);            
            hold on;
            grid on;
            set(gca,'FontSize',12,'fontWeight','bold');
            
            for objNum = 1:obj.objects 
                globalPosition = obj.position(:,objNum);
                globalVelocity = obj.velocity(:,objNum);
                globalQuaternion = obj.quaternion(:,objNum);
                % GET THE ROTATION MATRIX GOING FROM BODY-GLOBAL
                [R_GB,~] = obj.quaternionRotationMatrix(globalQuaternion); 
                
                % PLOT THE POSITIONS
                scatter3(globalPosition(1),globalPosition(2),globalPosition(3),'r');
                % PLOT THE VELOCITY VECTORS
                quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
                               globalVelocity(1),globalVelocity(2),globalVelocity(3),'c');
                quiv.AutoScaleFactor = 1;
                % ADD ANNOTATION
                annotationText = sprintf('   object %s',num2str(objNum));
                text(globalPosition(1),globalPosition(2),globalPosition(3),annotationText);
                
                % PLOT THE ROTATION TRIAD              
                % DEFINE TRIAD
                colours = {'r','g','b'};
                % REPLOT TRIAD
                for j = 1:size(obj.globalTriad,2)
                    % GET THE ROTATED TRIAD AXES
                    localVector = R_GB*obj.globalTriad(:,j);
                    % DRAW THE LOCAL TRIAD
                    quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
                                      localVector(1),localVector(2),localVector(3),colours{j}); 
                    quiv.AutoScaleFactor = 1;
                end
            end
            % Increment the plot counter
            plotnum = plotnum + 1; 
            % SAVE THE FIGURE TO THE WORKING DIRECTORY
            saveas(figureHandle,'scenario.png');
        end 
    end
    methods (Static)
        % PLOT SCENARIO FROM INITIALISED OBJECT-INDEX
        function [figureHandle] = plotObjectIndex(objectIndex)
            % This function is designed to plot the configured objectIndex 
            % using their global properties and thier simulation 
            % object-types
            
            global plotnum
            
            % DETERMINE PLOT PROPERTIES
            if ~exist('plotnum','var') || isempty(plotnum)
                plotnum = 1;    % Default to first plot
            end
                        
            % GENERATE THE FIGURE
            figureHandle = figure(plotnum);
            axis('equal');
            xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
            hold on; grid on;
            ax = gca;
            set(ax,'FontSize',12,'fontWeight','bold');
                       
            % OBJECT COUNTERS
            agents = 0; obstacles = 0; waypoints = 0;
            for index = 1:numel(objectIndex)
                % GET THE OBSTACLES GLOBAL PARAMETERS
                objectName       = objectIndex{index}.name;
                objectRadius     = objectIndex{index}.VIRTUAL.radius;
                objectType       = objectIndex{index}.VIRTUAL.type;
                globalPosition   = objectIndex{index}.VIRTUAL.globalPosition;
                globalVelocity   = objectIndex{index}.VIRTUAL.globalVelocity;
                globalQuaternion = objectIndex{index}.VIRTUAL.quaternion;
                % GET THE ROTATION MATRIX GOING FROM BODY-GLOBAL
                [R_q] = OMAS_geometry.quaternionToRotationMatrix(globalQuaternion);
                % PLOT THE OBJECT ORIENTATION TRIAD
                OMAS_graphics.drawTriad(figureHandle,globalPosition,R_q');
                % PLOT THE POSITIONS
                scatter3(globalPosition(1),globalPosition(2),globalPosition(3),'r');               
                % PLOT THE VELOCITY VECTORS
                quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
                               globalVelocity(1),globalVelocity(2),globalVelocity(3),'c');
                quiv.AutoScaleFactor = 1;
                % DETERMINE REPRESENTATION IF THE OBJECT HAS GEOMETRY
                if size(objectIndex{index}.GEOMETRY.vertices,1) < 1
                    % REPRESENT AS SPHERE WITH DEFINED RADIUS
                    [geometry] = OMAS_graphics.defineSphere(globalPosition,objectRadius);
                    vertexData = geometry.vertices;
                    faceData = geometry.faces;
                    entityFaceAlpha = 0.3;
                    entityLineWidth = 0.1;
                    entityEdgeAlpha = 0.1;                                 % Show representative shapes with higher alpha              
                else
                    vertexData = objectIndex{index}.GEOMETRY.vertices*R_q + globalPosition';
                    faceData = objectIndex{index}.GEOMETRY.faces;
                    entityFaceAlpha = 0.8;
                    entityLineWidth = 1;
                    entityEdgeAlpha = 0.8; %0.2; 
                end
                % REPRESENT GEOMETRY AS A PATCH
                entityHandle = patch(ax,...
                                    'Vertices',vertexData,...
                                    'Faces',faceData,...
                                    'EdgeColor','k',...
                                    'EdgeAlpha',entityEdgeAlpha,...
                                    'FaceAlpha',entityFaceAlpha,...
                                    'FaceLighting','gouraud',...
                                    'LineWidth',entityLineWidth);
                % PLOT REPRESENTATION
                switch objectType
                    case OMAS_objectType.agent
                        set(entityHandle,'FaceColor','b');
                        agents = agents + 1;
                    case OMAS_objectType.obstacle
                        set(entityHandle,'FaceColor','r');
                        obstacles = obstacles + 1;
                    case OMAS_objectType.waypoint
                        set(entityHandle,'FaceColor','g');
                        waypoints = waypoints + 1;
                    otherwise
                        set(entityHandle,'FaceColor','m');
                end
                % ADD ANNOTATION
                annotationText = sprintf('    %s [ID:%s]',objectName,num2str(objectIndex{index}.objectID));
                text(globalPosition(1),globalPosition(2),globalPosition(3),char(annotationText));
            end
            % ADD TITLE
            titleStr = sprintf('Test scenario: %.0f agents, %.0f obstacles and %.0f waypoints.',agents,obstacles,waypoints);
            title(titleStr);
            hold off;
        end
        % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
        function [config] = configurationParser(defaultConfig,scenarioParameters)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called from a get_scenario file.
            
            % INPUT HANDLING
            assert(mod(numel(scenarioParameters),2) == 0,' Please provide list of parameter:value pairs');
            % ASSIGN CONFIG TEMPLATE 
            config = defaultConfig;
            pairNum = numel(scenarioParameters)/2;
            for n = 1:pairNum
                % PULL PARAMETER/VALUE PAIRS
                parameterLabel = scenarioParameters{2*n - 1}; 
                parameterValue = scenarioParameters{2*n};
                % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                if isfield(config,parameterLabel)
                    config.(parameterLabel) = parameterValue;   % Make a substitution
                end
            end
        end
    end
end

