%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_3DVO.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_3DVO_qimeng < agent
    
    %% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % Tuning parameters
        obstacleSafetyFactor;               % Modify the apparent size of the obstacle
        % gif toggle
        giffOn = 0;
    end
    %%  CLASS METHODS
    methods
        %% CONSTRUCTOR
        function obj = agent_3DVO_qimeng(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
                namestr = '';                              % Assign default naming scheme
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(namestr);                           % Get super class 'agent'
            % AGENT SPECIFIC PARAMETERS
            obj.sensorRange = 50;                         % Virtual sensor range
            obj.obstacleSafetyFactor = 1.0;               % Modify the apparent size of the obstacle
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
            visualiseAgent = 1;
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
            
            % 2. BUILD AVOIDANCE PLANE
            for item = 1:length(obstacleSet)
                [AVO,DVO,alphaVO] = obj.defineVelocityObstacle(obstacleSet(item),visualiseProblem);
                VO(item) = struct('apex',AVO,'axis',DVO,'angle',alphaVO);
                % Elevation angle
                theta = obstacleSet(item).inclinationAngle;
                % Azimuth angle
                psi = obstacleSet(item).azimuthAngle;
                % the rotation angle of the plane (around x-axis)
                % choose a vertical plane in this case
                phi = 90*pi/180;
                % the free parameter of VO effective circular base
                beta = 0*pi/180;
                % the numbers of beta
                n = 40;
                % rotation matrice
                r_theta_psi = [cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0;0 0 1]*[cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
                r_phi = [1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
                for i = 1:n
                    % the free parameter of the generating lines
                    a = (VO(item).apex(3)*cos(phi) - VO(item).apex(2)*sin(phi))/((sin(phi)*cos(theta) + cos(phi))*sin(psi) + ((cos(phi)*cos(theta) - sin(phi)*sin(psi)*sin(theta))*sin(beta) - sin(phi)*cos(psi)*cos(beta))*tan(VO(item).angle));
                    % coordinates in collision cone
                    cc = r_theta_psi*[a;a*tan(VO(item).angle)*cos(beta);a*tan(VO(item).angle)*sin(beta)];
                    % coordinates in velocity obstacle
                    vo = cc + VO(item).apex;
                    % coordinates in avoidance plane
                    % a vertical plane in this case (phi = 90 degrees)
                    vo_phi = r_phi*vo;
                    % form a matrix of candidate velocites using coordinates
                    % of vo_phi
                    candidate_v(1,i) = vo_phi(1);
                    % For vertical plane, the second coordinate should be
                    % zero
                    candidate_v(2,i) = 0;
                    candidate_v(3,i) = vo_phi(3);
                    % calculate with a new beta next time
                    beta = beta + 2*pi/n;
                end
                % double the number of candidates as the plane is symmetric
                candidate_v(1:3,n+1:2*n) = -candidate_v(1:3,1:n);
            end
            
            % 3. CALCULATE THE ESCAPE VELOCITIES FROM THE CANDIDATES
            escapeVelocities = getEscapeVelocities(obj,candidate_v,VO);
            if isempty(escapeVelocities)
                warning('Agent %s has no viable escape velocities.',obj.name);
                escapeVelocities = [0;0;0];
            end
            % 4. PREPARE THE SEARCH SPACE
            searchMatrix = vertcat(escapeVelocities,zeros(1,size(escapeVelocities,2)));
            for i = 1:size(searchMatrix,2)
                searchMatrix(4,i) = sum(searchMatrix(1:3,i).^2);           % Add a dimension with the absolute vector value
            end
            
            % 5. DETERMINE THE MINIMUM (MAGNETUDE) ESCAPE VELOCITY AVAILABLE
            [~,I] = min(searchMatrix(4,:),[],2);                           % Return the index of the smallest vector
            
            % 6. USING THE SELECTED ESCAPE VELOCITY---
            %    CONVERT TO DYNAMIC ACCELERATION REQUEST
            desiredAcceleration = searchMatrix(1:3,I)/dt;
            
            % 7. PLOT THE AVOIDANCE SITUATION
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
            % 8. GET THE NEW STATE VECTOR
            newState = obj.stateDynamics(dt,desiredAcceleration,zeros(3,1));
            % 9. UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
    end
    
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        %% DRAW THE VELOCITY OBSTACLE PROBLEMS
        function [AVO,DVO,alphaVO] = defineVelocityObstacle(obj,obstacle,plotOn)
            % This function assembles the geometric collision avoidance
            % problem from the relative information observed from the
            % virtual sensor systems.
            
            %% DECLARE STATIC OBSTACLE KNOWLEDGE
            % Obstacle information
            Doi      = obstacle.position;                               % Relative position
            Vr   = obstacle.velocity;                               % Relative velocity
            mDoi  = sqrt(sum(Doi.^2));                            % Scalar relative seperation
            UDoi = obj.unit(Doi);                                % Unit relative position
            
            %% DEFINE THE ANGLE REFERENCE AXES
            unit_referenceAxis  = [0;0;1];                                 % Unit reference axis
            referenceVector = [1;0;0];                                     % The reference axis defaults to the agents direction
            problemAxis = cross(referenceVector,UDoi);
            if problemAxis == 0
                referenceVector = unit_referenceAxis;                       % Else form a comparison to the unit y axis
                problemAxis = cross(referenceVector,UDoi);           % Axes then align with the aerospace convention
            end
            
            %% DETERMINE THE VELOCITY OBSTACLE PROJECTION
            origin = [0;0;0];
            Dvo = origin + Doi + Vr;                          % Determine center of the projected sphere
            
            % radius of the protect zone
            rpz = obj.VIRTUAL.size + obstacle.size*obj.obstacleSafetyFactor;
            
            Rvo = -rpz*UDoi;                   % Unit radius of the obstacle (direction of doi)
            alphaVO = asin(rpz/mDoi);                  % Get rotation angle for tangential vector
            
            % DEFINE THE TANGENT VECTORS
            angle = (pi/2) - alphaVO;                                      % Define the leading tangent angle
            [tangRadiusL] = obj.rotateVectorAboutAxis(Rvo,problemAxis,angle);
            LT = Doi + tangRadiusL;                         % The leading tangent vector
            LTP = Dvo + tangRadiusL;                        % The leading tangent point
            
            %% DEFINE FINAL VELOCITY OBSTACLE ATTRIBUTES
            AVO = Vr;                                          % Define the VO origin point
            DVO = (dot(LT,Doi)/mDoi^2)*Doi;    % Define the tangent projection on AB
            
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
                [Xb,Yb,Zb] = obj.getSphere(Doi,obstacle.size);
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                    'FaceColor','g',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % DRAW THE VELOCITY OBSTACLE PROJECTION
                [X_vo,Y_vo,Z_vo] = obj.getSphere(Dvo,rpz);
                sphB = mesh(X_vo,Y_vo,Z_vo);
                set(sphB,'facealpha',0.2,...
                    'FaceColor','r',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % PLOT THE RELATIVE VELOCITY OF THE OBSTACLE
                q = quiver3(Doi(1),Doi(2),Doi(3),...
                    Vr(1), Vr(2), Vr(3),...   % Plot the obstacles velocity from its relative position
                    'b','filled');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENTIAL RADIUS OF THE VO
                q = quiver3(Dvo(1),Dvo(2),Dvo(3),...
                    tangRadiusL(1),tangRadiusL(2),tangRadiusL(3),'k');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(AVO(1),AVO(2),AVO(3),...
                    LT(1),LT(2),LT(3),'k');
                q.AutoScaleFactor = 1;
                % ASSEMBLE VELOCITY OBSTACLE CONES
                nodes = 20;
                [Cone] = obj.vectorCone(AVO,Dvo,LTP,nodes);
            end
        end
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
                    [flag] = obj.isInCone(VOlist(VOnumber).apex,unit_coneAxis,...
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
            Doi = pointB - pointA;
            mDoi = sqrt(sum((Doi).^2)); % Axis vector properties
            tangent = radialPoint - pointA;
            mod_tangent = sqrt(sum((tangent).^2));       % Tangental vector properties
            
            % DEFINE THE TANGENT-AXIS PROJECTION
            trueAB = (dot(tangent,Doi)/mDoi^2)*Doi;
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
            angle_X1X2 = acos(dot(unit_Vx,Doi)/(norm(unit_Vx)*mDoi))*180/pi;
            
            % Get rotation axis
            axis_rot = cross([1 0 0],Doi);
            
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
    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Static)
        % COMPARE POINT TO CONE GEOMETRY
        function [flag] = isInCone(apex,unit_coneAxis,mod_coneAxis,coneAngle,point)
            % INPUTS:
            % Avo
            % mod_VOaxis
            % unit_VOaxis
            % alpha      - Cone angle
            % probePoint - Test point
            
            % Compare the angle made between the probeVector and the
            % unitVO axis. Is less than alpha? else not in cone.
            
            probeVector = point - apex;           % Define validation vector
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
        % DRAW SPHERE
        function [X,Y,Z] = getSphere(position,radius)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius
            [X,Y,Z] = sphere(40);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
        end
        
        % ROTATE VECTOR THROUGH AN ANGLE, AROUND A GIVEN AXIS
        function [newVector] = rotateVectorAboutAxis(oldVector,Doi,theta)
            % This function is designed to calculate a vector
            % following a rotation around a given axis vector, through a
            % given angle.
            % INPUTS:
            % oldVector  - The initial vector
            % Doi - The axis of rotation
            % theta      - The angle of rotation
            % OUTPUTS:
            % newVector  - The rotated 3D vector
            
            % NORMALISE THE AXIS VECTOR
            Doi = Doi/(sqrt(sum(Doi.^2)));  % Normalize rotation axis
            
            % GET THE CROSS PRODUCT PROECTION BETWEEN THE AXIS AND VECTOR
            crossVector = cross(Doi,oldVector);
            
            % DETERMINE THE MAPPING OF EACH VECTOR COMPONENTS
            newVector = cos(theta)*oldVector ...
                + (crossVector)*sin(theta)  ...
                + Doi*(dot(Doi,oldVector))*(1 - cos(theta));
        end
        % UNIT VECTOR
        function [unitVector] = unit(inputVector)
            % This function returns the unit vector of a vector
            unitVector = inputVector/sqrt(sum(inputVector.^2));
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]