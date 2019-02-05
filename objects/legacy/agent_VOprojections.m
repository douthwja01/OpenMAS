%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_VOprojecions.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_VOprojections < agent
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties    
        % Performance Parameters
        maxAccelerations;
        minAccelerations;                   % Primative acceleration limits
        feasabliltyMatrix;                  % The stored achievable velocity matrix
        % Tuning parameters
        obstacleSafetyFactor;               % Modify the apparent size of the obstacle
        pointDensity;                       % Must be odd to allow a zero value
        
        % gif toggle
        giffOn = 0;
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR 
        function obj = agent_VOprojections(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
               namestr = '';                              % Assign default naming scheme
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(namestr);                           % Get super class 'agent'
            % AGENT SPECIFIC PARAMETERS
            obj.maxAccelerations = [300;300;3];
            obj.minAccelerations = [-300;-300;3];         % Primative acceleration limits
%             obj.maxAccelerations = [25;25;100];
%             obj.minAccelerations = [-25;-25;-100];      % Primative acceleration limits
            obj.sampleFrequency = 19;                     % Virtual computation frequency
            obj.sensorRange = 50;                         % Virtual sensor range
            obj.obstacleSafetyFactor = 1.0;               % Modify the apparent size of the obstacle
            obj.pointDensity = 7;                         % Must be odd to allow a zero value
            % INITIALISE THE FEASIBILIY MATRIX
            [obj.feasabliltyMatrix,~] = obj.getFeasabilityGrid(obj.minAccelerations,obj.maxAccelerations,obj.pointDensity);
        
            % VIRTUAL DEFINITION
            obj.VIRTUAL.detectionRange = obj.sensorRange; % Assign the range attribute to the SIM VIRTUAL attribute
        end
        %% AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
                                          
            % CALCULATE THE CONTROL LOOP CYCLE
            if 0 == rem(TIME.currentTime,(1/obj.sampleFrequency))               % WAIT FOR LOOP CYCLE SCHEDULE
                % GET THE NEW STATE VECTOR
                newState = obj.stateDynamics(dt,[0;0;0],[0;0;0]);
                % UPDATE THE CLASS GLOBAL PROPERTIES
                obj = obj.updateGlobalProperties(dt,newState);
                return
            end
            
            % CHECK FOR ENV UPDATE
            obstacleSet = varargin{1}; % The detected objects    
            
            % IF NO OBJECTS IN RANGE YET
            if isempty(obstacleSet)
                % GET THE NEW STATE VECTOR
                newState = obj.stateDynamics(dt,[0;0;0],[0;0;0]);
                % UPDATE THE CLASS GLOBAL PROPERTIES
                obj = obj.updateGlobalProperties(dt,newState);
                return
            end
            
            % PLOT THE LOCAL AVOIDANCE PROBLEM
            visualiseProblem = 0;
            visualiseAgent = 4;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis vis3d;
                view([-78 50]);
                xlabel('x_{m}');
                ylabel('y_{m}');
                zlabel('z_{m}');
                visualiseProblem = 1;
            end
            
            % 2. BUILD VELOCITY OBSTACLE SET
            for item = 1:length(obstacleSet)
                [VOorigin,VOaxis,VOangle] = obj.defineVelocityObstacle(obstacleSet(item),visualiseProblem);
                VO(item) = struct('origin',VOorigin,'axis',VOaxis,'angle',VOangle);
            end
            
            % 3. DEFINE REGION OF FEASBILITY
            availableVelocities = obj.feasabliltyMatrix*dt;                % Convert to velocity space            
            
            % 4. CALCULATE THE ESCAPE VELOCITIES FROM THE AVAILABLE
            escapeVelocities = getEscapeVelocities(obj,availableVelocities,VO);
            if isempty(escapeVelocities)
               warning('Agent %s has no viable escape velocities.',obj.name);
               escapeVelocities = [0;0;0];
            end
            % 5 PREPARE THE SEARCH SPACE
            searchMatrix = vertcat(escapeVelocities,zeros(1,size(escapeVelocities,2)));
            for i = 1:size(searchMatrix,2)
                searchMatrix(4,i) = sum(searchMatrix(1:3,i).^2);           % Add a dimension with the absolute vector value
            end
            
            % 6. DETERMINE THE MINIMUM (MAGNETUDE) ESCAPE VELOCITY AVAILABLE
            [~,I] = min(searchMatrix(4,:),[],2);                           % Return the index of the smallest vector          
            
            % 7. USING THE SELECTED ESCAPE VELOCITY---
            %    CONVERT TO DYNAMIC ACCELERATION REQUEST
            desiredAcceleration = searchMatrix(1:3,I)/dt;
            
            % 8. PLOT THE AVOIDANCE SITUATION
            if visualiseProblem && obj.objectID == visualiseAgent
                % PLOT THE ESCAPE VELOCITY PROBLEM
                scatter3(desiredAcceleration(1)*dt,...
                         desiredAcceleration(2)*dt,...
                         desiredAcceleration(3)*dt,'r','filled');
                %  PLOT THE LOCAL AGENT VELOCITY
                q = quiver3(0,0,0,desiredAcceleration(1)*dt...
                            ,desiredAcceleration(2)*dt...
                            ,desiredAcceleration(3)*dt,'k');
                q.AutoScaleFactor = 1;
%                 view([-90 0 ]);
                view([-45 45 ]);
                drawnow;
                % ANNOTATION WITH SIMTIME
                annoString = sprintf('Time: %ss',num2str(TIME.currentTime));
                dim = [0.85 0.1 0.1 0.1];
                annotation('textbox',dim,'String',annoString,'FitBoxToText','off');

                % DRAW AND COLLATE GIF FRAME
                im = frame2im(getframe(overHandle));
                [imind,cm] = rgb2ind(im,256);
                
                % APPEND THE FRAMES
                slomoFactor = 2;
                fname = sprintf('coneTransformation[%s]',obj.name);
                if ~obj.giffOn
                    imwrite(imind,cm,strcat(fname,'.gif')...
                        ,'gif', 'Loopcount',inf,'DelayTime',dt*slomoFactor);
                    obj.giffOn = 1;
                else
                    imwrite(imind,cm,strcat(fname,'.gif')...
                        ,'gif','WriteMode','append','DelayTime',dt*slomoFactor); % 0.05
                end
                hold off;
                close;
            end
            % 9. GET THE NEW STATE VECTOR
            newState = obj.stateDynamics(dt,desiredAcceleration,zeros(3,1));
            % 10. UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
    end
    
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        %% VELOCITY OBSTACLES METHOD
        function [escapeVelocities] = getEscapeVelocities(obj,potentialVelocities,VOlist)
            % This function takes a matrix of potential velocities and
            % compares them against the agents Velocity Obstacle set. A
            % matrix of escape velocities is then produced.
            % INPUTS:
            % potentialVelocities - The feasible velocity matrix
            % VOlist              - The velocity obstacle list
            % OUTPUTS:
            % escapeVelocities    - The available escape velocity matrix
            
            escapeVelocities = [];   % Prepare output container
            
            % FOR EACH VIABLE VELOCITY, CHECK AGAINST CONES
            for i = 1:size(potentialVelocities,2)
                % Nominate a 3D velocity coordinate
                candidatePoint = potentialVelocities(:,i); 
                % Check points against the VO set
                VOnumber = 1; flag = 0;
                while flag ~= 1 && VOnumber <= size(VOlist,2)   
                    % Calculate axis properties
                    mod_coneAxis = sqrt(sum(VOlist(VOnumber).axis.^2));
                    unit_coneAxis = VOlist(VOnumber).axis/mod_coneAxis;
                    
                    % DETERMINE WHETHER POINT IS WITHIN VO
                    [flag] = obj.isInCone(VOlist(VOnumber).origin,unit_coneAxis,...
                                          mod_coneAxis,VOlist(VOnumber).angle,...
                                          candidatePoint);
                    % Move to the next VO
                    VOnumber = VOnumber + 1;
                end
                % If the point lies outside all obstacles
                if ~flag 
                    escapeVelocities = horzcat(escapeVelocities,candidatePoint);
                end
            end
        end
        % DRAW THE VELOCITY OBSTACLE PROBLEMS
        function [VOorigin,VOaxis,VOangle] = defineVelocityObstacle(obj,obstacle,plotOn)
            % This function assembles the geometric collision avoidance
            % problem from the relative information observed from the
            % virtual sensor systems.
            
            %% DECLARE STATIC OBSTACLE KNOWLEDGE
            origin      = [0;0;0];   
            lambAB      = obstacle.observedPosition;                       % Relative position
            velocityB   = obstacle.observedVelocity;                       % Relative velocity
            mod_lambAB  = sqrt(sum(lambAB.^2));                            % Scalar relative seperation
            unit_lambAB = obj.unit(lambAB);                                % Unit relative position
            unit_referenceAxis  = [0;0;1];                                 % Unit reference axis
            
            %% DEFINE THE ANGLE REFERENCE AXES
%             referenceVector = unit_velocityB;                              % The reference axis defaults to the obstacles direction
            referenceVector = [1;0;0];  
            problemAxis = cross(referenceVector,unit_lambAB);
            if problemAxis == 0
               referenceVector = unit_referenceAxis;                       % Else form a comparison to the unit y axis
               problemAxis = cross(referenceVector,unit_lambAB);           % Axes then align with the aerospace convention
            end     
                 
            %% DETERMINE THE VELOCITY OBSTACLE PROJECTION
            configurationRadius = obj.VIRTUAL.size + obstacle.size*obj.obstacleSafetyFactor;
            % The projection is defined as the obstacle (in the
            % configuration space) augmented by its relative velocity.
            radiusVO = -configurationRadius*unit_lambAB;                   % Unit radius of the obstacle (direction of lambdaAB)
            VOangle = asin(configurationRadius/mod_lambAB);                  % Get rotation angle for tangential vector          
            
            % VELOCITY OBSTACLE VECTOR DESCRIPTION
            lambVO = origin + lambAB + velocityB;                          % Determine center of the projected sphere
                        
            % DEFINE THE TANGENT VECTORS
            angle = (pi/2) - VOangle;                                      % Define the leading tangent angle 
            [tangRadiusL] = rotateVectorAboutAxis(obj,radiusVO,problemAxis,angle);
            leadingTangent = lambAB + tangRadiusL;                         % The leading tangent vector
            leadingTangentPoint = lambVO + tangRadiusL;

            %% DEFINE FINAL VELOCITY OBSTACLE ATTRIBUTES
            VOorigin = velocityB;                                          % Define the VO origin point
            axisVector = lambAB;                                           % Define the VO axis vector
            mod_axisVector = mod_lambAB;                                   % Define the VO axis magnitude 
            VOaxis = (dot(leadingTangent,axisVector)/mod_axisVector^2)*axisVector; % Define the tangent projection on AB

            % PROBLEM VISUALISATION
            if plotOn
                % DRAW AGENT FORWARD VELOCITY
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(0,0,0,1,0,0,'r');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,-1,0,'g');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,0,-1,'b');
                q.AutoScaleFactor = 1;
                
                % DRAW OBSTACLE REPRESENTATION
                [Xb,Yb,Zb] = obj.getSphere(lambAB,obstacle.size);                            
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                         'FaceColor','g',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0); 
                % DRAW THE VELOCITY OBSTACLE PROJECTION
                [X_vo,Y_vo,Z_vo] = obj.getSphere(lambVO,configurationRadius);  
                sphB = mesh(X_vo,Y_vo,Z_vo); 
                set(sphB,'facealpha',0.2,...
                         'FaceColor','r',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0);
                % PLOT THE RELATIVE VELOCITY OF THE OBSTACLE
                q = quiver3(lambAB(1),lambAB(2),lambAB(3),...
                             velocityB(1), velocityB(2), velocityB(3),...   % Plot the obstacles velocity from its relative position
                            'b','filled');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENTIAL RADIUS OF THE VO
                q = quiver3(lambVO(1),lambVO(2),lambVO(3),...
                            tangRadiusL(1),tangRadiusL(2),tangRadiusL(3),'k');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(VOorigin(1),VOorigin(2),VOorigin(3),...
                            leadingTangent(1),leadingTangent(2),leadingTangent(3),'k');
                q.AutoScaleFactor = 1;
                % ASSEMBLE VELOCITY OBSTACLE CONES
                nodes = 10;
                [Cone] = obj.vectorCone(VOorigin,lambVO,leadingTangentPoint,nodes);
            end
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
            coneColour = 'r';
            coneEdgeColour = 'r';
            coneAlpha = 0.1;
            coneClosed = 0;
            coneLines = 0;
            
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
                     'FaceAlpha',coneAlpha,...
                     'EdgeColor',coneEdgeColour);        % Cone verticies
            if coneClosed==1
                set([EndPlate1 EndPlate2],...
                    'AmbientStrength',1,...
                    'FaceColor',coneColour,...
                    'FaceLighting','gouraud',...
                    'FaceAlpha',coneAlpha);         % End-plate 
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
            if obj.validateVector(oldVector,3,1) ~= obj.validateVector(axisVector,3,1)
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
    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Static)
        % COMPARE POINT TO CONE GEOMETRY
        function [flag] = isInCone(origin,unit_coneAxis,mod_coneAxis,coneAngle,point)      
            % INPUTS:
            % VOorigin
            % mod_VOaxis
            % unit_VOaxis
            % alpha      - Cone angle
            % probePoint - Test point 
            
            % Compare the angle made between the probeVector and the
            % unitVO axis. Is less than alpha? else not in cone.
            
            probeVector = point - origin;           % Define validation vector
            [unit_probeVector] = probeVector/sqrt(sum(probeVector.^2));
            theta = acos(dot(unit_probeVector,unit_coneAxis));      % The angle between cone axis and probe vector
            
            % CHECK POINT RELATION TO CONE GEOMETRY
            flag = 0;
            if theta <= coneAngle
                probeVOProjection = dot(probeVector,unit_coneAxis); %Calculate the projection of the probeVector onto VOaxis
                if probeVOProjection <= mod_coneAxis
                    % Is the projection of the point less than the cone length
                    % YES, IN CONE
                    flag = 1;
                end
            end
        end
        % BUILD A UNIT, CUBOID POINT CLOUD
        function [cubePoints,plotMesh] = getFeasabilityGrid(minVector,maxVector,pointDensity)
            % This function assembles a unit cuboid particle cloud of
            % potential velocity points. 
            % INPUTS:
            % minVector    - A vector of minimum accelerations
            % maxVector    - A vector of maximum accelerations
            % pointDensity - The number of descrete steps 
            % OUTPUTS:
            % cubePoints   - A vector of the 3D points within the cloud
            % plotMesh     - A mesh structure representing the volume
            
            % Get the descrete dimensional steps
            xDimension = linspace(minVector(1),maxVector(1),pointDensity);
            yDimension = linspace(minVector(2),maxVector(2),pointDensity);
            zDimension = linspace(minVector(3),maxVector(3),pointDensity);

            % Generate the 3D point list
            cubePoints = zeros(3,size(xDimension,2)^3);
            counter = 1;
            for i = 1:size(xDimension,2)
                for j = 1:size(xDimension,2)
                    for k = 1:size(xDimension,2)
                        cubePoints(1,counter) = xDimension(i);
                        cubePoints(2,counter) = yDimension(j);
                        cubePoints(3,counter) = zDimension(k);
                        counter = counter + 1;
                    end
                end
            end
            
            % GET THE MATLAB PLOT MESH
            getMesh = 0;
            if getMesh
                [plotMesh.xx,plotMesh.yy,plotMesh.zz] = meshgrid(xDimension,yDimension,zDimension);
            else 
                plotMesh = 0;
            end            
        end
        % DRAW SPHERE
        function [X,Y,Z] = getSphere(position,radius)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius
            [X,Y,Z] = sphere(40);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
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
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]