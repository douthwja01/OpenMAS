%% BASIC COLLISION AVOIDANCE AGENT (agent_CA.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed as a child of the agent class "agent", used to 
% define the properties of a generic Collision Avoidance agent.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef agent_CA0 < agent
%%% ARdrone CHILD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the ARdrone specific properties for importing a
    % small quadrotor into the environment.
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties
        % Sensor Parameters
        sensorRange = 80;
        
        % Performance Parameters
        maxAccelerations = [10;10;10];
        minAccelerations = [-10;-10;-10];
        
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR 
        function obj = agent_CA0(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
               namestr = ''; % Assign values to asset_args
            end
            % INITIALISE THE AGENT SUPERCLASS
            obj@agent(namestr);
            % DECLARE ANY INITIAL PARAMETERS
            obj.VIRTUAL.detectionRange = obj.sensorRange; % Range of 200m
            % VIRTUAL DEFINITION
        end
        
%%      UPDATE MECHANISMS
        % AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if ~isnumeric(varargin{1})
               error('Agent timestep is not numeric.');
            else
                dt = varargin{1};
            end
            
            % 1. OBSTACLE KNOWLEDGE
            obstacleSet = varargin{2}; % The detectable objects
            if isempty(obstacleSet) || length(varargin) < 2
                % NO DETECTIONS, CONTINUE ON CURRENT COURSE
                obj = updateStateVector(obj,[0;0;0],[0;0;0],dt);
                return
            end
            
            % 2. BUILD VELOCITY OBSTACLE SET
%             VOset = zeros(size(obstacleSet));
%              for item = 1:length(obstacleSet)
%                 [VObstacle] = obj.drawVelocityObstacle(obstacleSet{item});
%                 [VObstacle] = obj.defineVelocityObstacle(obstacleSet{item});
%                 VOset(item) = VObstacle;
%              end
            
            % 3. DETERMINE FEASIBLE ACCELERATIONS
%             [linearAcceleration,angularAcceleration] = achievableAccelerations();
%             maxVelocities = obj.maxAccelerations*dt;
%             minVelocities = obj.minAccelerations*dt;
            
            %[FRmatrix,FRregion] = obj.defineFeasibleRegion(maxVelocities,minVelocities);
            
            %trisurf(FRregion,FRmatrix(1,:)',FRmatrix(2,:)',FRmatrix(3,:)','Facecolor','red','FaceAlpha',0.1)
            %drawnow;
%             pause;
            % Need to add the current velocity to the change in velocity
            %pause;
            % ALLOCATE GLOBAL ACCELERATIONS
            linearAcceleration = [0;0;0];
            angularAcceleration = [0;0;0];
            
            % 4. UPDATE THE STATE VECTOR
            obj = updateStateVector(obj,linearAcceleration,angularAcceleration,dt);
        end
    end
    
%%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        %% VELOCITY OBSTACLES METHOD
        % VELOCITY OBSTACLE CONSTRUCTION
        function [VObstacle] = defineVelocityObstacle(obj,obstacle)
            % This function derives a volicty obstacle definition for an
            % agent within the detectable region. This algorithm selects
            % the position and velocity information from the detected.
            
            fprintf('|| Detection Object: %s, Obstacle: %s.\n',obj.name,obstacle.name);
                     
            % OBSTACLE OBSERVATION
            radiusB = obstacle.size + obj.VIRTUAL.size;
            
            % RELEVANT OBJECT INFORMATION
            relativeObstacleState = obstacle.state;            
            lambAB    = relativeObstacleState(1:3,1);   % Relative position
            velocityB = relativeObstacleState(4:6,1);   % Relative velocity
            unit_velocityB = velocityB/abs(sqrt(sum(velocityB.^2)));
            
            % VECTOR DEFINITIONS
            unit_velocityA = [1;0;0];
            mod_lambAB = abs(sqrt(sum(lambAB.^2)));     % Scalar relative seperation
            unit_lambAB = lambAB/mod_lambAB;            % Unit relative position
            VOorigin = velocityB;                       % Define the origin of the Velocity Obstacle 
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% DRAW THE VECTOR PROBLEM
            % The "reference" parameters refer the definition of the
            % problem plane. The problem is described in the plane of the
            % referenceAxis, which is defined by the arguement between the
            % seperation lambAB and a reference vector.
            
            referenceVector = unit_velocityB;
            referenceAxis = cross(referenceVector,unit_lambAB);
            if referenceAxis == 0
               referenceVector = unit_velocityA;                   % Else form a comparison to the unit y axis
               referenceAxis = cross(referenceVector,unit_lambAB); % Axes then align with the aerospace convention                
            end
            
            % Plot rotation axis acquire perpendicular radius vector
            rad_lambAB = -radiusB*unit_lambAB;  % Unit radius vector B in direction of A
            alpha = asin(radiusB/mod_lambAB);   % Get rotation angle for tangential vector
            
            % Plot the leading tangental radius
            angle = (pi/2) - alpha;
            [tangRadiusL] = rotateVectorAboutAxis(obj,rad_lambAB,referenceAxis,angle);
            
            % Plot the leading tangent line
            lambL = lambAB + tangRadiusL;
            
            %% GENERATE THE VELOCITY OBSTACLE CONE
            % The cone needs to be constructed from the positional vectors
            % along the new VO lambAB
            % VOorigin - The start of the cone
            % lambVO   - The center of th eprojected sphere
            originPoint = VOorigin;
            radialPoint = originPoint + lambL;
            axispoint = originPoint + lambAB;
            nodes = 10;
            [Cone] = obj.vectorCone(originPoint,axispoint,radialPoint,nodes);
            
            % Return the cone aas the velocity obstacle
            VObstacle = Cone;
        end
        % DRAW THE VELOCITY OBSTACLE PROBLEMS
        function [VObstacle] = drawVelocityObstacle(obj,obstacle)
            % This function derives a volicty obstacle definition for an
            % agent within the detectable region. This algorithm selects
            % the position and velocity information from the detected.
            
            fprintf('|| Detection Object: %s, Obstacle: %s.\n',obj.name,obstacle.name);
                     
            % OBSTACLE OBSERVATION
            radiusA = obj.VIRTUAL.size - obstacle.size; % Modify two radi
            radiusB = obstacle.size + obj.VIRTUAL.size;
            
            % RELEVANT OBJECT INFORMATION
            relativeObstacleState = obstacle.relativeState;            
            lambAB    = relativeObstacleState(1:3,1);   % Relative position
            velocityB = relativeObstacleState(4:6,1);   % Relative velocity
            unit_velocityB = velocityB/abs(sqrt(sum(velocityB.^2)));
            
            % VECTOR DEFINITIONS
            unit_xAxis = [1;0;0];
            velocityA = obj.velocity;
            mod_lambAB = abs(sqrt(sum(lambAB.^2)));     % Scalar relative seperation
            unit_lambAB = lambAB/mod_lambAB;            % Unit relative position
            VOorigin = velocityB;                       % Define the origin of the Velocity Obstacle 
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% PLOT THE SCENARIO
            % BUILD REFERENCE OBJECTS
            [X,Y,Z] = sphere(40);
            Xa = X; Ya = Y; Za = Z;
            Xb = X.*radiusB + lambAB(1);
            Yb = Y.*radiusB + lambAB(2);
            Zb = Z.*radiusB + lambAB(3);
            
            % PLOT CONFIGURATION
            f = figure(1);
            hold on;
            grid on;
            axis equal;
            xlabel('x_{m}');
            ylabel('y_{m}');
            zlabel('z_{m}');
            
            % GENERATE THE OBJECTS
            
            %sphA = mesh(Xa,Ya,Za); % CA Agent
            %set(sphA,'facealpha',0);
            sphB = mesh(Xb,Yb,Zb); % Obstacle
            set(sphB,'facealpha',0);
            
%%          PLOT PROBLEM VECTORS
            % Plot the seperation
%             quiver3(0,0,0,lambAB(1),lambAB(2),lambAB(3),'r','filled');
            % DRAW THE OBJECT VELOCITIES
            q = quiver3(0,0,0,unit_xAxis(1),unit_xAxis(2),unit_xAxis(3),'r','filled','LineWidth',3);
            q.AutoScaleFactor = 1;
            q = quiver3(lambAB(1),lambAB(2),lambAB(3),velocityB(1),velocityB(2),velocityB(3),'b','filled');
            q.AutoScaleFactor = 1;
            % Plot the origin of the velocity obstacle
            q = quiver3(0,0,0,VOorigin(1),VOorigin(2),VOorigin(3),'k','filled');
            q.AutoScaleFactor = 1;
            % Plot the cone axis
            q = quiver3(VOorigin(1),VOorigin(2),VOorigin(3),lambAB(1),lambAB(2),lambAB(3),'r','filled');
            q.AutoScaleFactor = 1;
            
            
            % Determine center of the projected sphere
            lambVO = VOorigin + lambAB;
            
            %% DRAW THE VECTOR PROBLEM
            % The "reference" parameters refer the definition of the
            % problem plane. The problem is described in the plane of the
            % referenceAxis, which is defined by the arguement between the
            % seperation lambAB and a reference vector.
            
            referenceVector = unit_velocityB;
            referenceAxis = cross(referenceVector,unit_lambAB);
            if referenceAxis == 0
               referenceVector = unit_xAxis;                   % Else form a comparison to the unit y axis
               referenceAxis = cross(referenceVector,unit_lambAB); % Axes then align with the aerospace convention
            end
            
            % Plot rotation axis acquire perpendicular radius vector
            quiver3(lambVO(1),lambVO(2),lambVO(3),referenceAxis(1),referenceAxis(2),referenceAxis(3),'b'); 
            
            rad_lambAB = -radiusB*unit_lambAB;  % Unit radius vector B in direction of A
            alpha = asin(radiusB/mod_lambAB);   % Get rotation angle for tangential vector
            
            % Plot the leading tangental radius
            angle = (pi/2) - alpha;
            [tangRadiusL] = rotateVectorAboutAxis(obj,rad_lambAB,referenceAxis,angle);
            q = quiver3(lambVO(1),lambVO(2),lambVO(3),tangRadiusL(1),tangRadiusL(2),tangRadiusL(3),'k');
            q.AutoScaleFactor = 1;
            
            % Plot the leading tangent line
            lambL = lambAB + tangRadiusL;
            q = quiver3(VOorigin(1),VOorigin(2),VOorigin(3),lambL(1),lambL(2),lambL(3),'c');
            q.AutoScaleFactor = 1;
            
            % Plot the trailing tangental radius
            angle = alpha - (pi/2);
            [tangRadiusT] = rotateVectorAboutAxis(obj,rad_lambAB,referenceAxis,angle);
            q = quiver3(lambVO(1),lambVO(2),lambVO(3),tangRadiusT(1),tangRadiusT(2),tangRadiusT(3),'k');
            q.AutoScaleFactor = 1;
            
            % Plot the trailing tangent
            lambT = lambAB + tangRadiusT;
            q = quiver3(VOorigin(1),VOorigin(2),VOorigin(3),lambT(1),lambT(2),lambT(3),'g');
            q.AutoScaleFactor = 1;
            
            %% GENERATE THE VELOCITY OBSTACLE CONE
            % The cone needs to be constructed from the positional vectors
            % along the new VO lambAB
            % VOorigin - The start of the cone
            % lambVO   - The center of th eprojected sphere
            originPoint = VOorigin;
            radialPoint = originPoint + lambT;
            axispoint = originPoint + lambAB;
            nodes = 10;
            [Cone] = obj.vectorCone(originPoint,axispoint,radialPoint,nodes);
            % Return the cone as the velocity obstacle
            VObstacle = Cone;
        end
        % DEFINE THE FEASIBLE REGION
        function [FRmatrix,FRregion] = defineFeasibleRegion(obj,maxVector,minVector)
           % This function is designed to define a scaled feasibility region.
           % INPUTS:
           % maxVector - The maximum accelerations at the next timestep
           % minVector - The minimum accelerations at the next timestep
           % OUTPUTS:
           % FRmatrix  - The vector definition of the cube
           
           % INPUT HANDLING
           if ~exist('minVector','var')
              minVector = [0;0;0]; 
           end
           
           % The unit vertices of the 3D feasible region 
           selectionMatrix = [0 1 0 0 1 1 0 1;...
                              0 0 1 0 0 1 1 1;...
                              0 0 0 1 1 0 1 1];
           % Create the relative maximum values
           maxVector = repmat((maxVector-minVector),1,8); 
           % Form the modified verticies of the feasible region
           FRmatrix = selectionMatrix.*maxVector + repmat(minVector,1,8);
           % GENERATE THE 3D REGION
           [FRregion,~] = boundary(FRmatrix(1,:)',FRmatrix(2,:)',FRmatrix(3,:)');
           
        end
        
        %% TOOLS & VECTOR OPERATIONS
        % CONSTRUCT CONE FROM VECTOR DEFINTIONS
        function [Cone] = vectorCone(obj,pointA,pointB,radialPoint,nodes)
            % This function is designed to construct a cone mesh between two 
            % 3D points Using a radial Point to specify the radial properties.
            % INPUTS:
            % pointA - First cone axis point (origin).
            % pointB - Second cone axis point.
            % radialPoint - Third point in space that defines the cones maximum radius
            % nodes  - The number of mesh nodes
            % OUTPUTS
            % Cone   - Cone mesh structure
            
            % INPUT HANDLING
            if ~exist('nodes','var')
                nodes = 10;
            end
            if ~exist('pointA','var')
                pointA = [0;0;0];
            end
            
            % MANUAL INPUTS
            coneColour = 'g';
            coneAlpha = 0.5;
            coneClosed = 1;
            coneLines = 1;
            
            % DEFINE VECTOR PROBLEM PARAMETERS
            axisVector = pointB - pointA;
            mod_axisVector = sqrt(sum((axisVector).^2)); % Axis vector properties
            tangent = radialPoint - pointA;
            mod_tangent = sqrt(sum((tangent).^2));       % Tangental vector properties
            
            % DEFINE THE TANGENT-AXIS PROJECTION
            trueAB = (dot(tangent,axisVector)/mod_axisVector^2)*axisVector;
            mod_trueAB = sqrt(sum((trueAB).^2));
            
            % GET THE RADIUS MODULUS
            mod_radius = sqrt(mod_tangent^2-mod_trueAB^2);
                                
            % Creating 2 circles in the YZ plane
            t=linspace(0,2*pi,nodes)';      % Create axis point set
            xa2 = zeros(length(t),1); 
            xa3 = zeros(size(xa2)); 
            xb2 = mod_radius*cos(t);
            xb3 = mod_radius*sin(t);        % Scale the second axis in the 
            
            % Creating the points in the X-Direction
            x1=[0 mod_trueAB];
            
            % Creating (Extruding) the cylinder points in the X-Directions
            xx1 = repmat(x1,length(xa2),1);
            xx2 = [xa2 xb2];  
            xx3 = [xa3 xb3]; % Concatinate second circle set
            
            % Drawing two filled cirlces to close the cylinder
            if coneClosed == 1
                EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
                EndPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'r');
            end
            
            % GENERATE THE CONE 
            % Plot the cone from the origin along x-axis and scale to size 
            Cone = mesh(real(xx1),real(xx2),real(xx3));
            
            % Get the planar rotation angle
            unit_Vx=[1 0 0];
            angle_X1X2 = acos(dot(unit_Vx,axisVector)/(norm(unit_Vx)*mod_axisVector))*180/pi;
            
            % Get rotation axis
            axis_rot = cross([1 0 0],axisVector); 
            
            % Rotating the plotted Cone and the end plate circles to the required
            % angles
            if angle_X1X2~=0 % Rotation is not needed if required direction is along X
                rotate(Cone,axis_rot,angle_X1X2,[0 0 0])
                if coneClosed==1
                    rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
                    rotate(EndPlate2,axis_rot,angle_X1X2,[0 0 0])
                end
            end
            
            % Till now Cone has only been aligned with the required direction, but
            % position starts from the origin. so it will now be shifted to the right
            % position
            if coneClosed == 1
                set(EndPlate1,'XData',get(EndPlate1,'XData') + pointA(1))
                set(EndPlate1,'YData',get(EndPlate1,'YData') + pointA(2))
                set(EndPlate1,'ZData',get(EndPlate1,'ZData') + pointA(3))
                
                set(EndPlate2,'XData',get(EndPlate2,'XData') + pointA(1))
                set(EndPlate2,'YData',get(EndPlate2,'YData') + pointA(2))
                set(EndPlate2,'ZData',get(EndPlate2,'ZData') + pointA(3))
            end
            set(Cone,'XData',get(Cone,'XData') + pointA(1))
            set(Cone,'YData',get(Cone,'YData') + pointA(2))
            set(Cone,'ZData',get(Cone,'ZData') + pointA(3))
            
            % SET THE COLOUR OF THE CONE AND END PLATES
            set(Cone,'AmbientStrength',1,...
                     'FaceColor',coneColour,...
                     'FaceLighting','gouraud',...
                     'FaceAlpha',coneAlpha);        % Cone verticies
            if coneClosed==1
                set([EndPlate1 EndPlate2],...
                    'AmbientStrength',1,...
                    'FaceColor',coneColour,...
                    'FaceLighting','gouraud',...
                    'FaceAlpha',coneAlpha);         % End-plate 
            else
                EndPlate1=[];
                EndPlate2=[];
            end
            
            % If lines are not needed making it disapear
            if coneLines == 0
                set(Cone,'EdgeAlpha',0)
            end
        end
        % ROTATE VECTOR THROUGH AN ANGLE, AROUND A GIVEN AXIS
        function [newVector] = rotateVectorAboutAxis(obj,oldVector,axisVector,theta)
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
            if obj.validateVector(oldVector,[3 1]) ~= obj.validateVector(axisVector,[3 1])
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
        % BUILD BODY FRAME REPRESENTATION
        function [angles,R]  = getVectorRotations(obj,referenceVector,inputVector)
            % This function is designed to calculate the angles between a
            % given reference vecto and an input vector and a rotation
            % matrix between them.
            
            rotationAxis = cross(referenceVector,inputVector); 
            rotationSkewMatrix = obj.getSkewSymmetricMatrix(rotationAxis);
            s = abs(sqrt(sum(rotationAxis.^2)));         % Sin of angle
            c = dot(referenceVector,inputVector);        % Cos of angle
            R = eye(3) + rotationSkewMatrix + rotationSkewMatrix^2*((1-c)/s^2);
            % Get angles
            angles = [0;asin(s);acos(c)];
            
        end
        % VECTOR ROTATION FUNCTION
        function [R,RX,RY,RZ] = getRotationMatrix(obj,phi,theta,psi)
           % This function is designed to rotate a given vector through the
           % angles given in the angle vector.
           
           % DECLARE THE INPUT ANGLES
%            phi   = angleVector(1); % rotation around the vector
%            theta = angleVector(2); % rotation in the xz plane (elevation)
%            psi   = angleVector(3); % rotation in the xy plane (azimuth)
           
           % DEFINE THE EULAR ROTATON MATRIX
           RX(1,1) = 1;
           RX(1,2) = 0;
           RX(1,3) = 0;
           RX(2,1) = 0;
           RX(2,2) = cos(phi);
           RX(2,3) = -sin(phi);
           RX(3,1) = 0;
           RX(3,2) = sin(phi);
           RX(3,3) = cos(phi);           
           
           RY(1,1) = cos(theta);
           RY(1,2) = 0;
           RY(1,3) = sin(theta);
           RY(2,1) = 0;
           RY(2,2) = 1;
           RY(2,3) = 0;
           RY(3,1) = -sin(theta);
           RY(3,2) = 0;
           RY(3,3) = cos(theta);
           
           RZ(1,1) = cos(psi);
           RZ(1,2) = -sin(psi);
           RZ(1,3) = 0;
           RZ(2,1) = sin(psi);
           RZ(2,2) = cos(psi);
           RZ(2,3) = 0;
           RZ(3,1) = 0;
           RZ(3,2) = 0;
           RZ(3,3) = 1;
           
           R = RZ*RY*RX;
           
           % DEFINE THE INVERSE ROTATION MATRIX
           Rinv = transpose(R);   % Invert for body to global
        end
        % CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
        function [outputMatrix] = getSkewSymmetricMatrix(obj,inputVector)
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
            
            outputMatrix = zeros(3,3);
            outputMatrix(1,2) = -inputVector(3);
            outputMatrix(1,3) =  inputVector(2);
            outputMatrix(2,1) =  inputVector(3);
            outputMatrix(2,3) = -inputVector(1);
            outputMatrix(3,1) = -inputVector(2);
            outputMatrix(3,2) =  inputVector(1); % Arrange the components
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]