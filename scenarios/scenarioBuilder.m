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
        quaternion;                             % The quaternion moving from the Global ENU to BODY ENU
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
        % PLANAR RING OF STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = planarRing(obj,varargin)
            % GENERATE ON CONSTANT ANGLE SCENARIO DEFINITION 
            % INPUTS:
            % pointA        - The first axis reference point.
            % pointB        - The second axis reference point; the ring center.
            % mod_radius    - The scalar radius of the ring
            % velMultiplier - The velocity multiplier for the object set.
            % OUTPUTS:
            % definition    - The scenario deifnition
            
            % DEFAULT CONDITIONS
            zeroAngle = 0;          % The rotation of the first node
            pointA = [0;0;-1];
            pointB = [0;0;0];
            radius = 10;
            velocityMultiplier = 1;
            
            % First of the axis points
            tmp = strncmpi(varargin,'pointA',6);
            if any(tmp)
                pointA = varargin{find(tmp) + 1};
            end
            % Second of the axis points
            tmp = strncmpi(varargin,'pointB',6);
            if any(tmp)
                pointB = varargin{find(tmp) + 1};
            end
            % The ring radius magnitude
            tmp = strncmpi(varargin,'radius',6);
            if any(tmp)
                radius = varargin{find(tmp) + 1};
            end
            % Angular position of the first node
            tmp = strncmpi(varargin,'zeroAngle',4);
            if any(tmp)
                zeroAngle = varargin{find(tmp) + 1};
            end
            % Velocity multiplier
            tmp = strncmpi(varargin,'velocityMultiplier',7); 
            if any(tmp)
                velocityMultiplier = varargin{find(tmp) + 1};
            end
            
            % DEFINE AN ANGLE FOR EQUAL DISTRIBUTION ABOUT THE RING
            nodalAngle = (2*pi)/obj.objects;          
            % HAND NODAL ANGLE TO ringAngle FUNCTION
            [ obj, scenarioConfig ] = obj.planarAngle('pointA',pointA,...
                                                      'pointB',pointB,...
                                                      'radius',radius,...
                                                      'zeroAngle',zeroAngle,...
                                                      'offsetAngle',nodalAngle,...
                                                      'velocityMultiplier',velocityMultiplier);
            % UPDATE SCENARIO
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
            
            % DEFAULT CONDITIONS
            zeroAngle = 0;          % The rotation of the first node
            pointA = [0;0;-1];
            pointB = [0;0;0];
            radius = 10;
            velocityMultiplier = 1;
            offsetAngle = pi/5;
                       
            % First of the axis points
            tmp = strncmpi(varargin,'pointA',6);
            if any(tmp)
                pointA = varargin{find(tmp) + 1};
            end
            % Second of the axis points
            tmp = strncmpi(varargin,'pointB',6);
            if any(tmp)
                pointB = varargin{find(tmp) + 1};
            end
            % The ring radius magnitude
            tmp = strncmpi(varargin,'radius',6);
            if any(tmp)
                radius = varargin{find(tmp) + 1};
            end
            % Angular position of the first node
            tmp = strncmpi(varargin,'zeroAngle',4);
            if any(tmp)
                zeroAngle = varargin{find(tmp) + 1};
            end
            % Nodal offset angle
            tmp = strncmpi(varargin,'offsetAngle',6);
            if any(tmp)
                offsetAngle = varargin{find(tmp) + 1};
            end
            % Velocity multiplier
            tmp = strncmpi(varargin,'velocityMultiplier',7); % |strncmpi(varargin,'velocityMultiplier',8); 
            if any(tmp)
                velocityMultiplier = varargin{find(tmp) + 1};
            end
            
            % UPDATE SCENARIO
            obj.name = sprintf('Planar offset angle of %srad.',num2str(offsetAngle));
            
            % GENERATE THE NODE(OBJECT) CARTAESIAN POSITIONS      
            % GET THE AXIS VECTOR PROPERTIES
            localXAxis = pointB - pointA;
            unit_axisVector = localXAxis/sqrt(sum(localXAxis.^2));
            
            % GET THE RADIAL VECTOR
            perpVector = cross(unit_axisVector,[1;0;0]);
            if sum(perpVector) == 0
                perpVector = cross(unit_axisVector,[0;1;0]);
            end
            perpVector = perpVector/sqrt(sum(perpVector.^2)); % Re-normalise the perpendicular vector
            % SCALE THE UNIT RADIAL VECTOR
            perpVector = radius*perpVector;               
            
            % GET THE NODE POINT SET
            nodalAngle = zeroAngle + pi;
            for node = 1:obj.objects
                % ROTATE THE RADIAL VECTOR TO DEFINE NODAL POSITIONS
                radialVector = obj.rotateVectorAboutAxis(perpVector,unit_axisVector,nodalAngle);
                % NODAL POSITIONS
                globalPosition = pointB + radialVector;
                % NODAL VELOCITIES
                unitVelocity = -(radialVector/sqrt(sum(radialVector.^2))); % -ve sign to make concentric velocities +ve             
                % NODAL ORIENTATIONS
                % We desire to orientate each object towards the center of
                % the ring. This is done by aligning the local body XY plane
                % With the plane of the ring. The local X axis is then
                % aligned with the cocentric global velocity. 
                                                
                % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                localTriad = zeros(3);
                localTriad(:,1) = obj.unit(unitVelocity);             % GET THE UNIT DIRECTION VECTOR
                localTriad(:,3) = obj.unit(unit_axisVector);          % UNIT RING AXIS
                localTriad(:,2) = obj.unit(cross(localTriad(:,1),localTriad(:,3)));
                
                % STORE THE OBJECTS IN SCENARIO DEFINITION 
                obj.position(:,node) = globalPosition;
                obj.velocity(:,node) = velocityMultiplier*unitVelocity;                
                obj.quaternion(:,node) = transpose(obj.getAnalyticalTriadRotation(obj.globalTriad,localTriad)); % GET THE GLOBAL-BODY QUATERNION
                % Object data is concatinated vertically
                % Increment the nodal position by the offset angle
                nodalAngle = nodalAngle - offsetAngle;
            end

            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);
            % UPDATE SCENARIO
            obj.name = 'Planar Angle.';
        end
        % SPHERICAL EQUAL DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = regularSphere(obj,varargin)
            % This function is designed to generate points equally
            % distributed about the circumference of a sphere.
            
            % DEFAULT CONDITIONS           
            sphereRadius = 10;
            sphereCenter = zeros(3,1);
            velocityMultiplier = 0;
            
            % POINT OF COCENTRICITY
            tmp = strncmpi(varargin,'center',3);
            if any(tmp)
                sphereCenter = varargin{find(tmp) + 1};
            end
            % SPHERE RADIUS
            tmp = strncmpi(varargin,'radius',3);
            if any(tmp)
                sphereRadius = varargin{find(tmp) + 1};
            end
            
            % VELOCITY MULTIPLIER
            tmp = strncmpi(varargin,'velocityMultiplier',7); 
            if any(tmp)
                velocityMultiplier = varargin{find(tmp) + 1};
            end
            
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
            XYZ = (pointCloud.*sphereRadius)';              
            XYZ = XYZ + sphereCenter;                                      % Offset and scale the sphere
            % SORT THE POINTS
            XYZ = unique(XYZ','rows');                                     % Remove repeat entries
            XYZ = sortrows(XYZ,1);
            XYZ = XYZ';
              
            planarIndices = any(XYZ == 0);
            nodalPositions = XYZ(:,planarIndices);
            nodalPositions = nodalPositions(:,1:obj.objects);
            
            % ALLOCATE RANDOM SET TO OBJECT GLOBAL POSITIONS
            for index = 1:obj.objects               
                % DETERMINE THE NODAL POSITIONS
                nodalPosition = sphereCenter + nodalPositions(:,index);
                
                errorFlag = 1;
                while errorFlag == 1   
                    % DEFINE GLOBAL VELOCITIES
                    unit_radii = obj.unit(sphereCenter - nodalPosition);
                    obj.velocity(:,index) = velocityMultiplier*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = obj.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = obj.unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.position(:,index) = nodalPosition;
                        % ROTATION
                        nodalQuaternion = obj.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                        obj.quaternion(:,index) = transpose(nodalQuaternion);
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

            obj.name = 'Equally spaced spherical distribution.';   
        end
        % COCENTRIC HELICAL STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = helix(obj,varargin)
            
            % DEFAULT CONFIGURATION
            defaultConfig = struct('file','scenario.mat',...
                                   'radius',10,...
                                   'waypointRadius',2,...
                                   'velocities',0,...
                                   'plot',0);
            % PARSE CONFIGURATIONS
            [config] = scenarioBuilder.configurationParser(defaultConfig,varargin)
            
            
            % SAVE THE RESULTING SCENARIO
            obj.name = 'Cocentric helix.';
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);   
        end
        
        % RANDOM SPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ obj,scenarioConfig ] = randomSphere(obj,varargin)
            % This fuction is designed to generate random positions 
            % distributed about the circumference of a sphere. 
            
            % DEFAULT CONDITIONS           
            sphereRadius = 10;
            sphereCenter = zeros(3,1);
            velocityMultiplier = 0;
            
            % POINT OF COCENTRICITY
            tmp = strncmpi(varargin,'center',3);
            if any(tmp)
                sphereCenter = varargin{find(tmp) + 1};
            end
            % SPHERE RADIUS
            tmp = strncmpi(varargin,'radius',3);
            if any(tmp)
                sphereRadius = varargin{find(tmp) + 1};
            end
            
            % VELOCITY MULTIPLIER
            tmp = strncmpi(varargin,'velocityMultiplier',7); 
            if any(tmp)
                velocityMultiplier = varargin{find(tmp) + 1};
            end
            
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
            XYZ = (pointCloud.*sphereRadius)';              
            XYZ = XYZ + sphereCenter;                                      % Offset and scale the sphere
            XYZ = unique(XYZ','rows');                                     % Remove repeat entries
            XYZ = XYZ';
            
            % RANDOMLY SAMPLE FROM THE POINT SET
            [nodalPositions,~] = datasample(XYZ,obj.objects,2,'Replace',false);   % Take one sample column from X coordinates

            % ALLOCATE RANDOM SET TO OBJECT GLOBAL POSITIONS
            for index = 1:obj.objects               
                % DETERMINE THE NODAL POSITIONS
                nodalPosition = sphereCenter + nodalPositions(:,index);
                
                errorFlag = 1;
                while errorFlag == 1   
                    % DEFINE GLOBAL VELOCITIES
                    unit_radii = obj.unit(sphereCenter - nodalPosition);
                    obj.velocity(:,index) = velocityMultiplier*unit_radii;     % Cocentric velocity
                    % GET THE LOCAL AXES (ROTATED IN THE GLOBAL FRAME)
                    localTriad = zeros(3);
                    localTriad(:,1) = unit_radii;                              % The body axis X direction
                    localTriad(:,3) = obj.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));      % The body axis Y direction
                    localTriad(:,2) = obj.unit(cross(localTriad(:,1),localTriad(:,3)));           % The body axis Z direction
                    % DETERMINE IF A QUATERNION ROTATION EXISTS
                    try
                        % ASSIGN POSITION
                        obj.position(:,index) = nodalPosition;
                        % ROTATION
                        nodalQuaternion = obj.getAnalyticalTriadRotation(obj.globalTriad,localTriad);
                        obj.quaternion(:,index) = transpose(nodalQuaternion);
                        errorFlag = 0;
                    catch
                        % IF ROTATION IS INVALID, RESAMPLE THE POSITION
                        [nodalPosition,~] = datasample(XYZ,1,2,'Replace',false);
                    end
                end
            end
            
            % SAVE THE RESULTING SCENARIO
            obj.name = 'Random cocentric sphere';
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);            
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
                localTriad(:,1) = obj.unit(normalObjectStates(velocityTuple,index));      % The body axis X direction
                localTriad(:,3) = obj.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = obj.unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternion(:,index) = transpose(obj.getAnalyticalTriadRotation(obj.globalTriad,localTriad));
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
                localTriad(:,1) = obj.unit(uniformObjectStates(velocityTuple,index));     % The body axis X direction
                localTriad(:,3) = obj.unit(cross(localTriad(:,1),obj.globalTriad(:,1)));  % The body axis Y direction
                localTriad(:,2) = obj.unit(cross(localTriad(:,1),localTriad(:,3)));       % The body axis Z direction
                obj.quaternion(:,index) = transpose(obj.getAnalyticalTriadRotation(obj.globalTriad,localTriad));
            end
            % SAVE THE RESULTING SCENARIO
            scenarioConfig = struct('objects',obj.objects,...
                                    'name',obj.name,...
                                    'position',obj.position,...
                                    'velocity',obj.velocity,...
                                    'quaternion',obj.quaternion);
%             save('scenario.mat','scenarioConfig');
        end
        
        %% SCENARIO PLOTTER
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
            
            %% GENERATE THE FIGURE
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
                
                %% PLOT THE ROTATION TRIAD              
                % DEFINE TRIAD
                colours = {'r','g','b'};
                %% REPLOT TRIAD
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
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods        
        % GET THE ROTATION ANGLES BETWEEN TWO VECTORS
        function [angles,R]  = getVectorRotations(obj,referenceVector,inputVector)
            % This function is designed to calculate the angles between a
            % given reference vector and an input vector; and a rotation
            % matrix between them.
            rotationAxis = cross(referenceVector,inputVector);
            rotationSkewMatrix = obj.getSkewMatrix(rotationAxis);
            s = abs(sqrt(sum(rotationAxis.^2)));         % Sin of angle
            c = dot(referenceVector,inputVector);        % Cos of angle
            R = eye(3) + rotationSkewMatrix + rotationSkewMatrix^2*((1-c)/s^2);
            % Get angles
            angles = abs([0;asin(s);acos(c)]);
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R_BG,R_GB] = quaternionRotationMatrix(obj,q)
            % This function defines the equivalent rotation matrix from the
            % qauternion q. Associated block "Quaternion to DCM"
            R_BG = zeros(3,3);           
            % Normalise the quaternion
            q = obj.unit(q);
            % Assemble the quaternion elements
            R_BG(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
            R_BG(1,2) = 2*(q(1)*q(4) + q(2)*q(3));
            R_BG(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
            R_BG(2,1) = 2*(q(3)*q(2) - q(1)*q(4));
            R_BG(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
            R_BG(2,3) = 2*(q(1)*q(2) + q(3)*q(4));
            R_BG(3,1) = 2*(q(1)*q(3) + q(2)*q(4));
            R_BG(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
            R_BG(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
            % Transpose the matrix for the global-body transformations
            R_GB = transpose(R_BG);
        end
        % ANALYTICALLY DERIVE THE QUATERNION FROM ONE TRIAD TO ANOTHER
        function [q] = getAnalyticalTriadRotation(obj,referenceTriad,targetTriad)
            % This function defines the quaternion describing the rotations
            % between a reference triad and a second triad.
            
            % NORMALISE THE TRIAD VECTORS
            for dim = 1:size(referenceTriad,2)
                [referenceTriad(:,dim)] = obj.unit(referenceTriad(:,dim));
                [targetTriad(:,dim)] = obj.unit(targetTriad(:,dim));
            end
            
            % EXTRACT REFERENCE AXES
            xAxis = targetTriad(:,1);
            zAxis = targetTriad(:,3);                                      % Get the rotated body axes
            xAxis_ref = referenceTriad(:,1);
            zAxis_ref = referenceTriad(:,3);                               % Get the reference unit ENU triad (trying to get to)
            
            % FIRST ALIGN THE Z-AXIS VECTORS //////////////////////////////
            % GET THE QUATERNION ROTATION TO ALLIGN THE Z-AXES
            [q_zAlign] = obj.vectorsToQuaternion(zAxis_ref,zAxis);         % Quaternion aligning global and body z vectors
            [R_zAlign] = obj.quaternionToRotationMatrix(q_zAlign);         % Equivalent rotation matrix
            % ALIGN THE X AXIS IN THE Z-PLANE
            xAxis_intermediate = R_zAlign*xAxis;
            % TAKE ITS PROJECTIONS IN THE XY PLANE & RENORMALISE
            xAxis_intermediate(3) = 0;
            [xAxis_intermediate] = obj.unit(xAxis_intermediate);
            % OTHERWISE JUST ALIGN THE X-AXIS VECTORS /////////////////////
            % GET THE QUATERNION ROTATION TO ALLIGN THE X-AXES
            [q_xAlign] = obj.vectorsToQuaternion(xAxis_ref,xAxis_intermediate);
            [R_xAlign,~] = obj.quaternionToRotationMatrix(q_xAlign);
            
            % COMPUTE THE COMPOSITE ROTATION MATRIX
            comp_rotation = R_xAlign * R_zAlign;
            % COVERT THE ROTATION MATRIX TO A QUATERNION DESCRIBING THE
            % ROTATION FROM TRIAD_REF TO TRIAD_FINAL
            q = rotm2quat(comp_rotation);
        end
        
        % OBJECT PLOTTING FUNCTIONS ///////////////////////////////////////
        % PLOT THE SCENARIO FROM THE OBJECT INDEX
        function [figureHandle] = plotObjectIndex(obj,objectIndex)
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
            set(gca,'FontSize',12,'fontWeight','bold');
            
            % THE GLOBAL TRIAD
%             globalTriad = [1,0,0;0,1,0;0,0,1];      % XYZ (Earth) frame of reference
            
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
                [R_q,~] = OMAS_axisTools.quaternionToRotationMatrix(globalQuaternion);
                % PLOT THE OBJECT ORIENTATION TRIAD
                colours = {'r','g','b'};
                % REPLOT TRIAD
                for j = 1:size(obj.globalTriad,2)
                    % GET THE ROTATED TRIAD AXES
                    localVector = R_q*obj.globalTriad(:,j);
                    % DRAW THE LOCAL TRIAD
                    quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
                        localVector(1),localVector(2),localVector(3),colours{j});
                    quiv.AutoScaleFactor = 1;
                end
                % PLOT THE POSITIONS
                scatter3(globalPosition(1),globalPosition(2),globalPosition(3),'r');
                % PLOT THE VELOCITY VECTORS
                quiv = quiver3(globalPosition(1),globalPosition(2),globalPosition(3),...
                               globalVelocity(1),globalVelocity(2),globalVelocity(3),'c');
                
                % ADD REPRESENTATIVE SPHERES
                switch objectType
                    case OMAS_objectType.agent
                        % IS AGENT
                        [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
                        objectColour = 'b';
                        agents = agents + 1;
                    case OMAS_objectType.obstacle
                        % IS OBSTACLE
                        [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
                        objectColour = 'r';
                        obstacles = obstacles + 1;
                    case OMAS_objectType.waypoint
                        % IS WAYPOINT
                        [X,Y,Z] = scenarioBuilder.drawSphere(globalPosition,objectRadius);
                        objectColour = 'g';
                        waypoints = waypoints + 1;
                    otherwise
                        error('[SCENARIO]\tObject type not recognised');
                end
                sphZone = mesh(X,Y,Z);
                set(sphZone,'facealpha',0.4,...
                    'FaceColor',objectColour,...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0.2);
             
                % ADD ANNOTATION
                annotationText = sprintf('    %s [ID:%s]',objectName,num2str(objectIndex{index}.objectID));
                text(globalPosition(1),globalPosition(2),globalPosition(3),char(annotationText));
            end
            % ADD TITLE
            titleStr = sprintf('Test scenario: %s agents, %s obstacles and %s waypoints.',num2str(agents),num2str(obstacles),num2str(waypoints));
            title(titleStr);
            hold off;
        end
    end
    
    methods (Static)
        %% GENERIC SCENARIO TOOLS /////////////////////////////////////////
        % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
        function [config] = configurationParser(defaultConfig,scenarioParameters)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called from a get_scenario file.
                       
            % MOVE THROUGH THE PARAMETER PAIRS ('agents',agentIndex)
            for parameterIndex = 1:numel(scenarioParameters)
                % FOR EACH USER-DEFINED PARAMETER
                givenParameter = scenarioParameters{parameterIndex};
                if ~ischar(givenParameter)
                    continue
                else
                    % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                    if isfield(defaultConfig,givenParameter)
                        defaultConfig.(givenParameter) = scenarioParameters{parameterIndex + 1};   % Make a substitution
                    end
                end
            end 
            % ERROR CHECKING //////////////////////////////////////////////
            % ADD SOME DEFAULT PARAMETERS TO THE 'defaultConfig'
            if ~isfield(defaultConfig,'file') || isempty(defaultConfig.file)
                defaultConfig.file = 'scenario.mat';
            end
            % CHECK AGENTS HAVE BEEN PROVIDED
            if ~isfield(defaultConfig,'agents') || isempty(defaultConfig.agents)
                error("[SCENARIO]\t You must specify and agent cell set using the 'agents' input parameter"); 
            end
            config = defaultConfig;
        end
        
        %% GENERIC PLOTTING FUNCTIONS /////////////////////////////////////

        % SPHERE DRAWING FUNCTION
        function [X,Y,Z] = drawSphere(position,radius)
            % BUILD REFERENCE SPHERE
            [X,Y,Z] = sphere(15);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
        end
        
        %% MATHMATICAL OPERATIONS /////////////////////////////////////////
        % GET THE QUATERNION BETWEEN TWO VECTORS
        function [q] = vectorsToQuaternion(u,v)
            q = zeros(4,1);
            % Normalise the quaternion
            u = u/norm(u);
            v = v/norm(v);
            % Get the axis vector
            q(2:4) = cross(u,v);
            % Define the rotation about that vector
            q(1) = sqrt((norm(u)^2)*(norm(v)^2)) + dot(u,v);
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R_AB,R_BA] = quaternionToRotationMatrix(q)
            % This function defines the equivalent rotation matrix from the
            % provided quaternion. (Associated block "Quaternion to DCM")
            % INPUT:
            % q    - The quaternion rotation
            % OUTPUT:
            % R_AB - The rotation matrix through A to B
            % R_BA - The reverse rotation matrix from B to A
            
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
            
            % INPUT HANDLING
            if size(oldVector,1) ~= 3 || size(axisVector,1) ~= 3
                display(oldVector); display(axisVector);
                warning('Input vector incorrectly defined.');
                newVector = [nan;nan;nan];
                return
            end
            
            % NORMALISE THE AXIS VECTOR
            axisVector = axisVector/(sqrt(sum(axisVector.^2)));  % Normalize rotation axis
            
            % GET THE CROSS PRODUCT PROECTION BETWEEN THE AXIS AND VECTOR
            crossVector = cross(axisVector,oldVector);
            
            % DETERMINE THE MAPPING OF EACH VECTOR COMPONENTS
            newVector = cos(theta)*oldVector ...
                + (crossVector)*sin(theta)  ...
                + axisVector*(dot(axisVector,oldVector))*(1 - cos(theta));
        end
        % CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
        function [outputMatrix] = getSkewMatrix(inputVector)
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
            unitVector = inputVector/norm(inputVector);
        end   
    end
end

