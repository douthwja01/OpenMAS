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
        name = 'Default Scenario.';
        % GLOBAL DATA
        % ALL GLOBAL DATA IS GIVEN IN THE EAST-NORTH-UP (ENU) FRAME OF
        % REFERENCE
        globalTriad = [1,0,0;0,1,0;0,0,1];      % XYZ (Earth) frame of reference
        positions;                              % Global positions
        velocities;                             % Global velocities
        quaternions;                            % The quaternion moving from the global XYZ
    end
    % PUBLIC METHODS
    methods  
        % CONSTRUCTION METHOD 
        function obj = scenarioBuilder(varargin)
            % The scenario builder is designed to define a scenario for a
            % given number of objects. An instance of the scenario builder
            % class must be created first to ensure the 
            
            % Access the OMAS-common properties
            addpath('environment/common');  % Access to common tools
            
            % Parse inputs against the builder
            [obj] = obj.configurationParser(obj,varargin);
   
            % Value checking
            assert(isnumeric(obj.objects),'Expecting a numeric value for field "objects"');
            assert(ischar(obj.name),'Expecting a string for field "name"');
            
            if obj.objects > 0
                % Initialise the scenario parameters
                obj.positions   = zeros(3,obj.objects);
                obj.velocities  = zeros(3,obj.objects);
                obj.quaternions = repmat([1;0;0;0],[1 obj.objects]);
            end
        end
               
        %% SCENARIO GENERATORS
        % PLANAR CONCENTRIC DISK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = planarDisk(obj,varargin)
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
            defaultConfig = struct('objects',0,...
                                   'pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'velocity',0,...
                                   'scale',1);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(isnumeric(config.pointA) && numel(config.pointA) == 3,'The reference point must be a 3D Cartesian point.'); 
            assert(isnumeric(config.pointB) && numel(config.pointB) == 3,'The centroid point must be a 3D Cartesian point.'); 
            assert(isnumeric(config.zeroAngle) && numel(config.zeroAngle) == 1,'The zero angle must be a scalar.'); 
            assert(isnumeric(config.velocity) && numel(config.velocity) == 1,'The object velocity must be a scalar.'); 
            assert(isnumeric(config.scale) && numel(config.scale) == 1,'The object padding scale must be a scalar.'); 
            
            % ////////////// DESIGN THE DISK DISTRIBUTION /////////////////
            % Update scenario label
            obj.name = 'Distributed planar disk.';
            % Define object number
            obj.objects = config.objects;
            % If there is only one object, simply place it at the center
            if obj.objects == 1
                obj.positions(:,1) = config.pointB;
                obj.velocities(:,1) = [1;0;0]*config.velocity;
                obj.quaternions(:,1) = [1;0;0;0];
                return
            end
            remainingObjects = obj.objects;
            layer = 1; 
            while layer < obj.objects && remainingObjects > 0
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
                [layerConfig] = obj.planarRing(...
                    'objects',objectsInLayer,...
                    'pointA',config.pointA,...
                    'pointB',config.pointB,...
                    'zeroAngle',config.zeroAngle,...
                    'radius',layerRadius);
                
                % EXTRACT THE LAYER GLOBAL PROPERTIES
                first = 1 + obj.objects - remainingObjects;
                last  = (first - 1) + objectsInLayer;
                indexVector = first:last;
                obj.positions(:,indexVector) = layerConfig.positions(:,1:objectsInLayer);
                obj.velocities(:,indexVector) = layerConfig.velocities(:,1:objectsInLayer);
                obj.quaternions(:,indexVector) = layerConfig.quaternions(:,1:objectsInLayer);
                % INCREMENT LAYER
                remainingObjects = remainingObjects - objectsInLayer;
                layer = layer + 1;
            end

        end
        % PLANAR RING OF STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = planarRing(obj,varargin)
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
            defaultConfig = struct('objects',0,...
                                   'pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'radius',10,...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'offsetAngle',pi/5,...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % Update scenario label
            obj.name = 'Planar ring.';
            % Define object number
            obj.objects = config.objects;
            
            % HAND NODAL ANGLE TO ringAngle FUNCTION
            [ obj ] = obj.planarAngle(...
                'objects',config.objects,...
                'pointA',config.pointA,...
                'pointB',config.pointB,...
                'radius',config.radius,...
                'zeroAngle',config.zeroAngle,...
                'offsetAngle',(2*pi)/config.objects,... % Define an angle for equal distribution
                'velocity',config.velocity);
        end
        % CONSTANT OFFSET ANGLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = planarAngle(obj,varargin)
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
            defaultConfig = struct('objects',0,...
                                   'pointA',[0;0;-1],...
                                   'pointB',[0;0;0],...
                                   'radius',10,...
                                   'zeroAngle',0,...                       % The rotation of the first node
                                   'offsetAngle',pi/5,...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(isnumeric(config.pointA) && numel(config.pointA) == 3,'The reference point must be a 3D Cartesian point.'); 
            assert(isnumeric(config.pointB) && numel(config.pointB) == 3,'The center point must be a 3D Cartesian point.'); 
            assert(isnumeric(config.radius) && numel(config.radius) == 1,'The arc radius must be a scalar.'); 
            assert(isnumeric(config.zeroAngle) && numel(config.zeroAngle) == 1,'The initial radial angle must be a scalar.'); 
            assert(isnumeric(config.offsetAngle) && numel(config.offsetAngle) == 1,'The offset radial angle must be a scalar.'); 
            assert(isnumeric(config.velocity) && numel(config.velocity) == 1,'The object velocity must be a scalar.');
            
            % Update scenario label
            obj.name = sprintf('Planar offset angle of %srad.',num2str(config.offsetAngle));
            % Define object number
            obj.objects = config.objects;
            
            % GENERATE THE NODE(OBJECT) CARTAESIAN POSITIONS      
            % GET THE AXIS VECTOR PROPERTIES
            localXAxis = config.pointB - config.pointA;
            unit_axisVector = unit(localXAxis);
            % GET THE RADIAL VECTOR
            perpVector = cross(unit_axisVector,[1;0;0]);
            if sum(perpVector) == 0
                perpVector = cross(unit_axisVector,[0;1;0]);
            end
            perpVector = unit(perpVector); % Re-normalise the perpendicular vector
            % SCALE THE UNIT RADIAL VECTOR
            perpVector = config.radius*perpVector;               
            
            % GET THE NODE POINT SET
            nodalAngle = config.zeroAngle + pi/2; % pi/2 aligns the first object with the x-axis
            for node = 1:config.objects
                % ROTATE THE RADIAL VECTOR TO DEFINE NODAL POSITIONS
                radialVector = OMAS_geometry.rodriguesRotation(perpVector,unit_axisVector,nodalAngle);
                % NODAL POSITIONS
                globalPosition = config.pointB + radialVector;
                % NODAL VELOCITIES
                unitVelocity = - unit(radialVector);        % -ve sign to make concentric velocities +ve             
                % NODAL ORIENTATIONS
                % We desire to orientate each object towards the center of
                % the ring. This is done by aligning the local body XY plane
                % With the plane of the ring. The local X axis is then
                % aligned with the cocentric global velocity. 
                                                
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                localTriad = zeros(3);
                localTriad(:,1) = unit(unitVelocity);             % GET THE UNIT DIRECTION VECTOR
                localTriad(:,3) = unit(unit_axisVector);          % UNIT RING AXIS
                localTriad(:,2) = unit(cross(localTriad(:,1),localTriad(:,3)));
                
                % STORE THE OBJECTS IN SCENARIO DEFINITION 
                obj.positions(:,node)   = globalPosition;
                obj.velocities(:,node)  = config.velocity*unitVelocity; 
                obj.quaternions(:,node) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                
                % Object data is concatinated vertically
                % Increment the nodal position by the offset angle
                nodalAngle = nodalAngle - config.offsetAngle;
            end
        end
        % SPHERICAL EQUAL DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = regularSphere(obj,varargin)
            % This function is designed to generate points equally
            % distributed about the circumference of a sphere.
                       
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',0,...
                                   'radius',10,...
                                   'center',[0;0;0],...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(isnumeric(config.radius) && numel(config.radius) == 1,'The sphere radius must be a scalar.'); 
            assert(isnumeric(config.center) && numel(config.center) == 3,'The sphere center must be a 3D Cartesian vector.'); 
            assert(isnumeric(config.velocity) && numel(config.velocity) == 1,'The object velocity must be a scalar.'); 
            
            % ATTEMPT TO ADD THE ICOSAHEDRALS LIBRARY
            try
                libraryPath = strcat(cd,'\scenarios\icosahedrals');
                addpath(libraryPath);
            catch
                error('Unable to load the icosahedral generator function, has it been moved?');
            end                 
            
            % UPDATE SCENARIO LABEL
            obj.name = 'Equally spaced spherical distribution.';  
            obj.objects = config.objects;
            
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
            nodalPositions = XYZ();
            nodalPositions = nodalPositions(:,1:obj.objects);
            % ALLOCATE RANDOM SET TO OBJECT GLOBAL POSITIONS
            for index = 1:obj.objects               
                % DETERMINE THE NODAL POSITIONS
                nodalPosition = config.center + nodalPositions(:,index);
                
                errorFlag = 1;
                while errorFlag == 1   
                    % DEFINE GLOBAL VELOCITIES
                    unit_radii = unit(config.center - nodalPosition);
                    obj.velocities(:,index) = config.velocity*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.positions(:,index) = nodalPosition;
                        % ROTATION
                        obj.quaternions(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                        errorFlag = 0;
                    catch
                        % IF ROTATION IS INVALID, RESAMPLE THE POSITION
                        [nodalPosition,~] = datasample(XYZ,1,2,'Replace',false);
                    end
                end
            end
        end
        % COCENTRIC HELICAL STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = helix(obj,varargin)
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',0,...
                                   'pointA',[0;0;0],...
                                   'pointB',[0;0;10],...
                                   'radius',10,...
                                   'offsetAngle',pi/5,...
                                   'velocities',0,...
                                   'plot',0);
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(ischar(config.file),'The provided file path must be a string.'); 
            assert(isnumeric(config.pointA) && numel(config.pointA) == 3,'The origin point must be a 3D Cartesian point.');
            assert(isnumeric(config.pointB) && numel(config.pointB) == 3,'The terminal point must be a 3D Cartesian point.');
            assert(isnumeric(config.radius) && numel(config.radius) == 1,'The helical radius must be a scalar.');
            assert(isnumeric(config.offsetAngle) && numel(config.offsetAngle) == 1,'The helical angle must be a scalar.');
            assert(isnumeric(config.velocity) && numel(config.velocity) == 1,'The object velocity must be a scalar.');
            
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
                obj.positions(:,i) = planarConfig.position(:,1);
                obj.velocities(:,i) = planarConfig.velocity(:,1);
                obj.quaternions(:,i) = planarConfig.quaternion(:,1);
                % SHIFT UP THE AXIS VECTOR
                pointA_i = pointB_i;    
            end
            
            obj.name = 'Cocentric helix.'; 
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.positions,...
                                    'velocity',obj.velocities,...
                                    'quaternion',obj.quaternions); 
        end
        % COLINEAR DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = line(obj,varargin)
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',0,...
                                   'pointA',[0;0;0],...
                                   'pointB',[0;0;10],...
                                   'heading',[1;0;0],...
                                   'velocities',1,...
                                   'plot',0);
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);

            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(isnumeric(config.pointA) && numel(config.pointA) == 3,'Please provide a 3D point to the start of the line.'); 
            assert(isnumeric(config.pointB) && numel(config.pointB) == 3,'Please provide a 3D point to the end of the line.');
            assert(isnumeric(config.heading) && numel(config.heading) == 3,'Please provide a 3D vector defining the object heading');
            assert(isnumeric(config.velocities) && numel(config.velocities) == 1,'Please provide a scalar initial speed along the heading vector');
            
            % Label the scenario
            obj.name = 'Apposing line configuration';            
            obj.objects = config.objects;

            % DISTRIBUTE THE OBJECTS EQUALLY ALONG THE AXIS
            if isnumeric(config.objects)
                nObjects = config.objects;
            else
                nObjects = numel(config.objects);
            end
            
            headingReference = unit(config.heading);    % Ensure the heading is a unit vector
            lineLength = norm(config.pointB - config.pointA);
            lineAxis = (config.pointB - config.pointA)/lineLength;
            % Linearly interpolate
            segments = linspace(0,lineLength,nObjects);
                        
            % Design the parameters for each object
            for node = 1:nObjects
                % OBJECT POSITION
                obj.positions(:,node)  = config.pointA + segments(node)*lineAxis;
                % OBJECT VELOCITY
                obj.velocities(:,node) = config.velocities*config.heading;
                % ORIENTATION
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)                  
                localTriad = zeros(3);
                localTriad(:,1) = headingReference;                        % GET THE UNIT DIRECTION VECTOR
                localTriad(:,3) = unit_mex(obj.globalTriad(:,3)); % GLOBAL Z-AXIS
                localTriad(:,2) = unit_mex(cross(localTriad(:,1),localTriad(:,3)));
                % RESOLVE ROTATION QUATERNION
                obj.quaternions(:,node) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end
        end
        
        % RANDOM SPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = randomSphere(obj,varargin)
            % This fuction is designed to generate random positions 
            % distributed about the circumference of a sphere. 
                        
            % DEFAULT CONFIGURATION
            defaultConfig = struct('objects',0,....
                                   'radius',10,...
                                   'center',[0;0;0],...
                                   'velocity',0);   
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % INPUT SANITY CHECK
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');
            assert(isnumeric(config.radius) && numel(config.radius) == 1,'Sphere radius must be a scalar.');
            assert(isnumeric(config.center) && numel(config.ceneter) == 3,'Sphere center must be a 3D Cartesian vector.');
            assert(isnumeric(config.velocity) && numel(config.velocity) == 1,'Agent velocity must be a scalar.');
            
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
                    unit_radii = unit(sphereCenter - nodalPosition);
                    obj.velocities(:,index) = config.velocity*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.positions(:,index) = nodalPosition;
                        % ROTATION
                        obj.quaternions(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
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
                                    'position',obj.positions,...
                                    'velocity',obj.velocities,...
                                    'quaternion',obj.quaternions);   
            % UPDATE SCENARIO LABEL
            obj.name = 'Random cocentric sphere';
        end  
        % DEFAULT RANDOM STATE GENERATOR - NORMAL %%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = random(obj,varargin)
            % This function is a simple placeholder to facilitate the use
            % of .random for quick access to normally distributed state sets.
            [ obj ] = obj.randomNormal(varargin);
        end
        % RANDOM (NORMAL DISTRIBUTION) SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = randomNormal(obj,varargin)
            % This function generates a random state vector set using a normal distribution. 
            % assuming the order of the states is [x;y;z;u;v;w;phi;theta;psi;p;q;r]
            % INPUT:
            % obj          - The scenario definition object
            % OUTPUT:
            % randomStates - The state vector set [states x objects]
            % obj          - The returned scenario object
                        
            % DEFAULT CONFIGURATION
            defaultConfig = struct(...
                'objects',0,....
                'velocity',0,...
                'positionGain',10,...
                'velocityGain',10,...
                'poseGain',pi);
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % Input sanity check
            assert(isnumeric(config.objects) && config.objects > 0,'Expecting a scalar number of objects greater than zero.');

            % DEFINE A RANDOM STATE PARAMETERS
            obj.name = 'Normal-randomised';
            % The number of objects
            obj.objects = config.objects;
            
            stateNum = 9;
            positionTuple = 1:3;
            velocityTuple = 4:6;
            angleTuple = 7:9;
            stateGains = [config.positionGain,config.velocityGain,config.poseGain];
            normalObjectStates = zeros(stateNum,1);
            
            for index = 1:obj.objects
                % GENERATE RANDOM NUMBER SET
                maxR = max(randn(stateNum,1));
                normalObjectStates(:,index) = randn(stateNum,1)/maxR; % returns an n-by-n matrix of random numbers.
                
                % ALLOCATE STATE VALUES
                normalObjectStates(positionTuple,index) = stateGains(1)*normalObjectStates(positionTuple,index);
                normalObjectStates(velocityTuple,index) = stateGains(2)*normalObjectStates(velocityTuple,index);
                
                % NODAL POSITIONS
                obj.positions(:,index)  = normalObjectStates(positionTuple,index);
                % NODAL VELOCITIES
                obj.velocities(:,index) = normalObjectStates(velocityTuple,index);
                                                    
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
                localTriad(:,1) = unit(normalObjectStates(velocityTuple,index));      % The body axis X direction
                localTriad(:,3) = unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternions(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end

        end
        % RANDOM (UNIFORM DISTRIBUTION) SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj ] = randomUniform(obj,varargin)
            % This function generates a random state vector set using a uniform distribution. 
            % assuming the order of the states is [x;y;z;u;v;w;phi;theta;psi;p;q;r]
            % INPUT:
            % obj          - The scenario definition object
            % OUTPUT:
            % obj          - The returned scenario object
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct(...
                'objects',0,....
                'velocity',0,...
                'positionGain',10,...
                'velocityGain',10,...
                'poseGain',pi);
            
            % PARSE CONFIGURATIONS
            [config] = obj.configurationParser(defaultConfig,varargin);
            
            % DEFINE A RANDOM STATE PARAMETERS
            obj.name = 'Uniformly-randomised';
            % The number of objects
            obj.objects = config.objects;
            
            stateNum = 9;
            positionTuple = 1:3;
            velocityTuple = 4:6;
            
            angleTuple = 7:9;
            stateGains = [config.positionGain,config.velocityGain,config.poseGain];
            uniformObjectStates = zeros(stateNum,1);
            
            for index = 1:obj.objects
                % GENERATE RANDOM NUMBER SET
                uniformObjectStates(:,index) = (rand([stateNum,1])-0.5).*2;       % returns an n-by-n matrix of random numbers.
                % ALLOCATE LINEAR STATE VALUES
                uniformObjectStates(positionTuple,index) = stateGains(1)*uniformObjectStates(positionTuple,index);
                uniformObjectStates(velocityTuple,index) = stateGains(2)*uniformObjectStates(velocityTuple,index);
                
                % NODAL POSITIONS
                obj.positions(:,index)  = uniformObjectStates(positionTuple,index);
                % NODAL VELOCITIES
                obj.velocities(:,index) = uniformObjectStates(velocityTuple,index);
                                
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
                localTriad(:,1) = unit(uniformObjectStates(velocityTuple,index));     % The body axis X direction
                localTriad(:,3) = unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternions(:,index) = OMAS_geometry.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
            end
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
            axis vis3d;
            view([45,45]);
            xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
            hold on; grid on;
            ax = gca;
            set(ax,'FontSize',12,'fontWeight','bold');
                       
            % OBJECT COUNTERS
            agents = 0; obstacles = 0; waypoints = 0;
            for index = 1:numel(objectIndex)
                % GET THE OBSTACLES GLOBAL PARAMETERS
                GLB = objectIndex{index}.GetGLOBAL();
                objectName  = objectIndex{index}.name;
                objectID    = objectIndex{index}.objectID;

                % GET THE ROTATION MATRIX GOING FROM BODY-GLOBAL
                [R_q] = OMAS_geometry.quaternionToRotationMatrix(GLB.quaternion);
                % PLOT THE OBJECT ORIENTATION TRIAD
                OMAS_graphics.drawTriad(figureHandle,GLB.position,R_q');
                % PLOT THE POSITIONS
                scatter3(GLB.position(1),GLB.position(2),GLB.position(3),'r');               
                % PLOT THE VELOCITY VECTORS
                quiv = quiver3(GLB.position(1),GLB.position(2),GLB.position(3),...
                               GLB.velocity(1),GLB.velocity(2),GLB.velocity(3),'c');
                quiv.AutoScaleFactor = 1;
                % DETERMINE REPRESENTATION IF THE OBJECT HAS GEOMETRY
                if size(objectIndex{index}.GEOMETRY.vertices,1) < 1
                    % REPRESENT AS SPHERE WITH DEFINED RADIUS
                    [geometry] = OMAS_graphics.defineSphere(GLB.position,GLB.radius);
                    vertexData = geometry.vertices;
                    faceData = geometry.faces;
                    entityFaceAlpha = 0.3;
                    entityLineWidth = 0.05;
                    entityEdgeAlpha = 0.1;                                 % Show representative shapes with higher alpha              
                else
                    vertexData = objectIndex{index}.GEOMETRY.vertices*R_q + GLB.position';
                    faceData   = objectIndex{index}.GEOMETRY.faces;
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
                switch GLB.type
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
                annotationText = sprintf('    %s [ID:%s]',objectName,num2str(objectID));
                text(GLB.position(1),GLB.position(2),GLB.position(3),char(annotationText));
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
            
            % Call the generic parameter overrider
            [config] = GetParameterOverrides_recursive(defaultConfig,scenarioParameters);
        end
    end
end

