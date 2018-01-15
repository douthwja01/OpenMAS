%% HYBRID RECIPROCAL VELOCITY OBSTACLE METHOD (agent_HRVO.m) %%%%%%%%%%%%%%
% This agent classe implements the Hybrid-Reciprocal Velocity Obstacle
% (HRVO) method outlined in: "The Hybrid Reciprocal Velocity Obstacle" by
% Jur van den Berg et al.

% Author: James A. Douthwaite

classdef agent_HRVO < agent_RVO
    
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE HRVO METHOD
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_HRVO(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent_RVO(varargin);                                       % Get the supercalss
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
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
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + obj.objectID);
                hold on; grid on;
                axis vis3d;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            observationSet = varargin{1};                                  % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(observationSet);
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredSpeed = 2;
            if ~isempty(obj.targetWaypoint)
                % DEFAULT THE HEADING VECTOR TOWARDS THE WAYPOINT
                waypointPosition = obj.targetWaypoint.state(1:3);
                desiredVelocity  = (waypointPosition/norm(waypointPosition))*desiredSpeed;
            else
                desiredVelocity = [1;0;0]*desiredSpeed; 
            end
                        
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % /////////////// PARSE CONTROLLER INPUTS /////////////////////
            targetSpeed = norm(desiredVelocity);
            targetHeading = desiredVelocity/targetSpeed;                   % Get the control inputs
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            useController = 0;
            if useController
                [d_heading,d_speed,obj] = obj.trajectoryController(targetHeading,targetSpeed);  
            else
                [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);
            end                     
            newSpeed = (norm(obj.localState(7:9)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = d_heading;                                        % Heading change relative to straight ahead
            
            % //////////////// APPLY KINEMATIC LIMITS /////////////////////
            [newHeading,newSpeed] = obj.kinematicContraints(dt,newHeading,newSpeed);
            newVelocity = [1;0;0]*newSpeed;                                % Convert to local velocity vector
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            [newState] = obj.stateDynamics_velocities(dt,newVelocity,newHeading);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                newState = obj.freezeAgent();
            end
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:length(obj.priorError),TIME.currentStep) = [newSpeed;newState(4:6)];         % Record the control inputs 
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
        end

        %% THE RECIPROCAL AVOIDANCE OVERRIDING FUNCTIONS
        % GET THE VO VELOCITY CORRECTION
        function [avoidance_V] = getAvoidanceCorrection(obj,dt,desired_V,knownObstacles,visualiseProblem)
            % This function calculates the collision avoidance heading and
            % velocity correction.
            
            % AGENT PROPERTIES
            current_V = obj.localState(7:9,1);
            capable_V = current_V + obj.feasabliltyMatrix*dt;              % Get all the feasible velocities this timestep 
                       
            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:length(knownObstacles)
                % DEFINE THE VELOCITY OBSTACLE PROPERTIES
                [VOorigin,VOaxis,VOangle] = obj.defineHybridReciprocalVelocityObstacle(knownObstacles(item),visualiseProblem);
                VO = vertcat(VO,struct('origin',VOorigin,'axis',VOaxis,'angle',VOangle));
            end
            
            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
            [escape_V] = obj.getEscapeVelocities(capable_V,VO);              % The viable velocities from the capable 
            if isempty(escape_V)
                warning('No viable escape velocities');
                avoidance_V = desired_V;
                return
            end
            
            % SEARCH THE VIABLE ESCAPE VELOCITIES FOR THE VELOCITY WITH THE 
            % SMALLEST DEVIATION FROM THE DESIRED VELOCITY          
            searchMatrix = vertcat(escape_V,zeros(1,size(escape_V,2)));
            for i = 1:size(searchMatrix,2)
                searchMatrix(4,i) = norm(desired_V - searchMatrix(1:3,i)); % Add a dimension with the absolute vector value
            end
            % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
            [~,minIndex] = min(searchMatrix(4,:),[],2);                    % Return the index of the smallest vector
            avoidance_V = searchMatrix(1:3,minIndex);
            
            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem
                scatter3(capable_V(1,:),capable_V(2,:),capable_V(3,:),'r'); 
                scatter3(escape_V(1,:),escape_V(2,:),escape_V(3,:),'g'); 
                q = quiver3(0,0,0,current_V(1),current_V(2),current_V(3),'b');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,desired_V(1),desired_V(2),desired_V(3),'g');
                q.AutoScaleFactor = 1;               
                q = quiver3(0,0,0,avoidance_V(1),avoidance_V(2),avoidance_V(3),'r');
                q.AutoScaleFactor = 1;
                scatter3(avoidance_V(1,:),avoidance_V(2,:),avoidance_V(3,:),'b','filled'); 
            end
        end
    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Access = public)
        % DEFINE THE RECIPCOCAL VELOCITY OBSTACLE
        function [HRVOorigin,HRVOaxis,HRVOgamma] = defineHybridReciprocalVelocityObstacle(obj,obstacle,plotOn)
            % This function designs a Hybrid-Reciprocal Velocity Obstacle
            % (HRVO) for a sensed obstacle. This method is that presented
            % in the paper "The Hybrid Reciprocal Velocity Obstacle" by
            % Jamie Snape et. al.
            % INPUTS:
            % obj      - The agent object
            % obstacle - The sensed obstacle from the agents memory
            % plotOn   - The output plot flag
            % OUTPUTS:
            % VOapex   - The apex of HRVO
            % VOaxis   - The vector describing the cone axis
            % gamma_VO - The cone angle from the axis vector
            
            % AGENT KNOWLEDGE
            v_a = obj.localState(7:9);
            r_a = obj.VIRTUAL.size;
            % OBSTACLE KNOWLEDGE
            origin    = [0;0;0];
            r_b       = obstacle.size;
            lambda_ab = obstacle.state(1:3);
            c_b       = obstacle.state(4:6);
            v_b = v_a + c_b;
            r_c = r_a + obj.obstacleSafetyFactor*r_b;
                        
            % ////////// PARAMETERS DEFINING THE PROBLEM PLANE ////////////
            mod_lambda_b  = sqrt(sum(lambda_ab.^2));                        % Scalar relative seperation
            unit_lambda_b = obj.unit(lambda_ab);                            % Unit relative position
            unit_referenceAxis  = [0;0;1];                                 % Unit reference axis
            
            % DEFINE THE PLANE OF AVOIDANCE                        
            referenceVector = [1;0;0];                                     % The reference axis defaults to the agents direction
            problemAxis = cross(referenceVector,unit_lambda_b);
            if problemAxis == 0
               referenceVector = unit_referenceAxis;                       % Else form a comparison to the unit y axis
               problemAxis = cross(referenceVector,unit_lambda_b);         % Axes then align with the aerospace convention
            end    
            
            % /////////// HYBRID-RECIPROCAL VELOCITY OBSTACLE /////////////
            % CONSTRUCT THE VELOCITY OBSTACLE
            r_VO = -r_c*unit_lambda_b;                                     % Unit radius of the obstacle (direction of lambdaAB)
            HRVOgamma = asin(r_c/mod_lambda_b);                            % Get rotation angle for tangential vector
            
            % ////// DEFINE THE LEADING & TRAILING) TANGENT VECTORS ///////
            leadingAngle = (pi/2) - HRVOgamma;                             % Define the angle to leading tangent point
            [leadingTangentialRadius]  = obj.rotateVectorAboutAxis(r_VO,problemAxis,leadingAngle);
            leadingTangentVector = lambda_ab + leadingTangentialRadius;    % The leading tangent vector
            
            trailingAngle = HRVOgamma - (pi/2);                             % Define the angle to leading tangent point
            [trailingTangentialRadius]  = obj.rotateVectorAboutAxis(r_VO,problemAxis,trailingAngle);
            trailingTangentVector = lambda_ab + trailingTangentialRadius;  % The leading tangent vector  
            
            % ///////// CALCULATE THE HRVO APEX/ORIGIN POSITION ///////////
            VOorigin = v_b;
            RVOorigin = (v_a + v_b)/2;
            % CALCULATE THE VECTOR INTERSECTIONS USING KNOW PARAMETERS
            PA = [RVOorigin,VOorigin]';                                    % Rays starting at the two origins
            VOtrailing = VOorigin + trailingTangentVector;              
            RVOleading = RVOorigin + leadingTangentVector;
            PB = [RVOleading,VOtrailing]';                                 % Ray termination point
            [P_intersect] = obj.intersect3Dvectors(PA,PB);
            HRVOoffset = P_intersect';

            % TRANSLATE THE COLLISION CONE INTO THE HRVO 
            HRVOorigin = origin + HRVOoffset;                              % Position of the HRVO origin
            HRVOaxis = (dot(leadingTangentVector,lambda_ab)/mod_lambda_b^2)*lambda_ab; % Define the tangent projection on AB
            HRVOlambda = HRVOorigin + lambda_ab;                           % The position of the projected sphere 
            HRVOtangentPoint = HRVOlambda + leadingTangentialRadius;       % The leading tangent of the HRVO    

            % PROBLEM VISUALISATION
            if plotOn
                hold on;
                % PLOTTING PROPERTIES               
                coneCenter = HRVOorigin + HRVOaxis;                           % Cone-tangent modified axis
                
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(0,0,0,1,0,0,'r');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,-1,0,'g');
                q.AutoScaleFactor = 1;
                q = quiver3(0,0,0,0,0,-1,'b');
                q.AutoScaleFactor = 1;
                
                % DRAW OBSTACLE REPRESENTATION
                [Xb,Yb,Zb] = obj.getSphere(lambda_ab,r_b);                            
                sphB = mesh(Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                         'FaceColor','k',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0); 
                % PLOT THE RELATIVE VELOCITY OF THE OBSTACLE
                q = quiver3(lambda_ab(1),lambda_ab(2),lambda_ab(3),...
                                 v_b(1),v_b(2),v_b(3),'b','filled');
                q.AutoScaleFactor = 1; 
                
                % \\\\\\  DRAW THE VELOCITY OBSTACLE PROJECTION \\\\\\\\\\\
                % PLOT THE VO ORIGIN
                q = quiver3(origin(1),origin(2),origin(3),...
                            HRVOorigin(1),HRVOorigin(2),HRVOorigin(3),'b');
                q.AutoScaleFactor = 1;
                % DRAW THE VELOCITY OBSTACLE PROJECTION
                [X_vo,Y_vo,Z_vo] = obj.getSphere(HRVOlambda,r_c);  
                sphB = mesh(X_vo,Y_vo,Z_vo); 
                set(sphB,'facealpha',0.2,...
                         'FaceColor','b',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0);
                % PLOT THE LEADING TANGENTIAL RADIUS
                q = quiver3(HRVOlambda(1),HRVOlambda(2),HRVOlambda(3),...
                            leadingTangentialRadius(1),leadingTangentialRadius(2),leadingTangentialRadius(3),'k');
                q.AutoScaleFactor = 1;
                % PLOT THE LEADING TANGENT VECTOR
                q = quiver3(HRVOorigin(1),HRVOorigin(2),HRVOorigin(3),...
                    leadingTangentVector(1),leadingTangentVector(2),leadingTangentVector(3),'k');
                q.AutoScaleFactor = 1;
                
                % ASSEMBLE VELOCITY OBSTACLE CONES
                nodes = 11;
                obj.vectorCone(HRVOorigin,coneCenter,HRVOtangentPoint,nodes,'b');
            end
            
        end
    end
    methods (Static)
        % INTERSECT 3D CARTESIAN VECTORS
        function [P_intersect] = intersect3Dvectors(PA,PB)
            % Find intersection point of lines in 3D space, in the least squares sense.
            % PA :          Nx3-matrix containing starting point of N lines
            % PB :          Nx3-matrix containing end point of N lines
            % P_Intersect : Best intersection point of the N lines, in least squares sense.
            % distances   : Distances from intersection point to the input lines
            % Anders Eikenes, 2012
            
            Si = PB - PA; %N lines described as vectors
            ni = Si ./ (sqrt(sum(Si.^2,2))*ones(1,3)); %Normalize vectors
            nx = ni(:,1); ny = ni(:,2); nz = ni(:,3);
            SXX = sum(nx.^2-1);
            SYY = sum(ny.^2-1);
            SZZ = sum(nz.^2-1);
            SXY = sum(nx.*ny);
            SXZ = sum(nx.*nz);
            SYZ = sum(ny.*nz);
            S = [SXX SXY SXZ;SXY SYY SYZ;SXZ SYZ SZZ];
            CX  = sum(PA(:,1).*(nx.^2-1) + PA(:,2).*(nx.*ny)  + PA(:,3).*(nx.*nz));
            CY  = sum(PA(:,1).*(nx.*ny)  + PA(:,2).*(ny.^2-1) + PA(:,3).*(ny.*nz));
            CZ  = sum(PA(:,1).*(nx.*nz)  + PA(:,2).*(ny.*nz)  + PA(:,3).*(nz.^2-1));
            C   = [CX;CY;CZ];
            P_intersect = (S\C)';
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]