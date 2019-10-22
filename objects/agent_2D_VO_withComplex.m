%% 2D GEOMETRIC COLLISION AVOIDANCE AGENT (agent_2D_VO.m) %%%%%%%%%%%%%%%%%

% Author: James A. Douthwaite

classdef agent_2D_VO_withComplex < agent_2D & agent_VO
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE 2D VO METHOD
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [obj] = agent_2D_VO_withComplex(varargin)
            % Construct the agent object and initialise with the following
            % specific parameters.
            
            % Call the super class
            obj@agent_2D(varargin); 

            % Omit superclass field
            obj.feasabliltyMatrix = [];                                    % Omit parent field
            
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin); 
            % Re-affirm associated properties   
            [obj] = obj.SetRadius(obj.radius);                             % Reaffirm radius against .VIRTUAL
            [obj] = obj.SetDetectionRadius(obj.neighbourDist);             % Define the sensor range by the neighbour distance
            % /////////////////////////////////////////////////////////////
        end
        % Setup - X = [x y psi dx dy dpsi]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            % INITIALISE THE 2D STATE VECTOR WITH CONCANTINATED VELOCITIES
            [obj] = obj.initialise_2DVelocities(localXYZVelocity,localXYZrotations);
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % obj      - The agent object
            % ENV      - The current environmental data
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated object
            
            % GET THE TIMESTEP
            dt = ENV.dt;
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet,waypointSet] = obj.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE

            % DEBUG PLOT
            visualiseProblem = 1;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(1);
                hold on; grid on;
                axis square vis3d;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
                view([0 90]);
                % GET THE OBJECT SCENE
                [overHandle] = obj.GetObjectScene(overHandle);
                [overHandle] = obj.GetAnimationFrame(ENV,overHandle,'objectScene.gif');
            end      
            
            % /////////// BEHAVIOUR OF THE AGENT_VO AGENT ////////////////
            
            % DEFAULT BEHAVIOUR 
            desiredSpeed = obj.nominalSpeed;
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            desiredHeadingVector = obj.GetTargetHeading();
            desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
                        
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet;agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = obj.GetAvoidanceCorrection(dt,desiredVelocity,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            
            % /////// COMPUTE STATE CHANGE FROM CONTROL INPUTS ////////////
            [obj] = obj.controller(dt,desiredVelocity);

            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'dx (m/s)','dy (m/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = obj.localState(4:6);         % Record the control inputs
        end
    end
    %% /////////////////////// AUXILLARY METHODS ///////////////////////////
    methods
        % CALCULATE THE NECESSARY 2D AVOIDANCE VELOCITY
        function [headingVector,speed] = GetAvoidanceCorrection(obj,dt,desiredVelocity,visualiseProblem)
            % This function calculates the 2D avoidance velocity command
            % and returns it to be achieved by the controller.
            
            % INPUT HANDLING
            assert(numel(desiredVelocity) == 2,'Desired velocity must be 2D');
            
            % AGENT KNOWLEDGE
            [p_i,v_i,r_i] = obj.GetAgentMeasurements(); % Its own position, velocity and radius
            
            % Define the obstacle list
            obstacleIDs = [obj.MEMORY([obj.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
            
            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:numel(obstacleIDs)
                % Get object data from memory structure
                p_j = obj.GetLastMeasurementByID(obstacleIDs(item),'position');
                v_j = obj.GetLastMeasurementByID(obstacleIDs(item),'velocity');
                
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;            % Maximum number of neighbours
%                 neighbourConditionB = norm(knownObstacles(item).position) < obj.neighbourDist;  % [CONFIRMED]
                neighbourConditionC = ~any(isnan(v_j));                    % Wait for a valid velocity reading
                % || ~neighbourConditionB
                if ~neighbourConditionA  || ~neighbourConditionC
                    continue
                end
                
                % OBSTACLE KNOWLEDGE
                r_j = obj.GetLastMeasurementByID(obstacleIDs(item),'radius');
                g_j = obj.GetLastMeasurementByID(obstacleIDs(item),'geometry');
                p_j = p_j + p_i;                                           % The relative -> absolute coordinates
                v_j = v_j(1:2,1) + v_i;                                    % Convert relative parameters to absolute
                tau_j = 0;
                
                % DEFINE THE VELOCITY OBSTACLE PROPERTIES
                if OMAS_objectType.agent == obj.GetLastMeasurementByID(obstacleIDs(item),'type')
                    % CALCULATE THE VELOCITY FOR COMPLEX VELOCITY
                    [VO_i] = obj.define2DVelocityObstacle_complex(p_i,v_i,r_i,p_j,v_j,g_j,tau_j);
                else
                    [VO_i] = obj.define2DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                end
                VO_i.objectID = obstacleIDs(item);                        % Add a unique identifier
                VO = [VO,VO_i];
            end
            
            % THE CLEAR PATH STRATEGY
            [avoidanceVelocity] = obj.strategy_clearPath(v_a,desiredVelocity,VO,visualiseProblem);
            
            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0];  % Retain previous heading
            end
        end
        % DECOMPOSE THE OBSTACLE GEOMETRIES INTO A VO SET
        function [VO] = define2DVelocityObstacle_complex(obj,p_a,v_a,r_a,p_b,v_b,geometry,tau)
            % This function computes the VO set from a list of complex
            % obstacles defined by their .GEOMETRY parameter.
            
            % NOTES:
            % - The .geometry parameter is the geometry of the object pre-
            %   positioned relative to the agent, rotated and scaled.
            
            fig = gcf;            
            hold on;
            axis equal;
            
            % INPUT HANDLING
            if ~obj.VIRTUAL.is3D
                p_a = [p_a;0]; v_a = [v_a;0];
                p_b = [p_b;0]; v_b = [v_b;0];
            end
            
            % CONSTANTS
            origin   = zeros(3,1);
            VOorigin = v_b;
            lambda   = p_b - p_a;
            unit_reference = [1;0;0];
            unit_normal = [0;0;1];
            
            % GET THE VERTICES PROJECTED ONTO THE X/Y PLANE
            [planarVertices] = OMAS_graphics.geometryPlanarProjection(unit_normal,geometry.centroid',geometry.vertices);
            
            % MOVE THROUGH THE PLANAR PROJECTIONS 
            vDot = zeros(size(planarVertices,1),1);
            onRightLogicals = zeros(size(planarVertices,1),1);
            for i = 1:size(planarVertices,1)
                % THE VEREX VECTOR
                vertexVector = planarVertices(i,:)' - origin;
                unit_vertex = vertexVector/norm(vertexVector);
                % ON RIGHT FLAG
                onRightLogicals(i) = -sign(dot(cross(unit_reference,unit_vertex),unit_normal));
                
                % THE VERTEX PROJECTION PERPENDICULAR TO THE REFERENCE 
                refProjection =  dot(unit_reference,vertexVector)*unit_reference;
                
                perpVector = vertexVector - refProjection
                
                vDot(i) = norm(perpVector);
                
%                 b = vertexVector + origin;
                q = quiver3(refProjection(1),...
                    refProjection(2),...
                    refProjection(3),...
                    perpVector(1),...
                    perpVector(2),...
                    perpVector(3),'-');
                q.AutoScaleFactor = 1;
            end
            
            
            
            % GET THE SIGNED LEFT/RIGHT PROJECTIONS 
            vertexProjections = vDot.*onRightLogicals;
            
            leftVertex = vertexProjections(onRightLogicals < 0)
            
            % LOGICALS
%             [~,I] = min(vDot,[],1)
            
            
            % THE DIRECTION OF ROTATION FROM THE FORWARD VECTOR
            % centroidIsRightLogical = -sign(dot(cross(unit_reference,unit_lambda_planar),planeNormal)); % 1 (right) -1 (left) % Direction of the rotation relative to the centerline
            
            % CENTER-LINE
%             q = quiver3(gca,origin(1),origin(2),origin(3),lambda_planar(1),lambda_planar(2),lambda_planar(3),'-.');
            
            % /////////////////////////////////////////////////////////////
%             for i = 1:size(lambda_vertices,1)
%                 % GET THE VERTEX PARAMETERS
%                 vertex_planar = lambda_vertices(i,:)';
%                 norm_vertex = norm(vertex_planar);
%                 unit_vertex = vertex_planar/norm_vertex;
%                 
%                 % VERTEX IS ON RIGHT/LEFT OF THE FORWARD VECTOR
%                 vertexIsOnRightLogical = -sign(dot(cross(unit_reference,unit_vertex),planeNormal));       % Direction of the rotation relative to the centerline
%                 
%                 % THE MINIMUM VECTOR FROM THE CURRENT VERTEX
%                 %                 perpVector = r_a*cross(planeNormal,unit_vertex);
%                 
%                 % /////////////////////////////////////////////////////////
%                 
%                 % THE VERTEX MODIFIER (r_a) ANGLE
%                 %                 alphaAngles(i) = asin(r_a/norm_vertex);
%                 vertexAngle = -vertexIsOnRightLogical*acos(dot(unit_reference,unit_vertex));
%                 if abs(vertexAngle) == pi
%                     %                     vertexAngle = 0;
%                     continue
%                 end
%                 % BUILD VERTEX ANALYSIS STRUCTURE
%                 vertexData(i) = struct('point',vertex_planar,...
%                                         'distance',norm_vertex,...
%                                         'perp',r_a*cross(planeNormal,unit_vertex),...
%                                         'isOnRight',vertexIsOnRightLogical,...
%                                         'angle',vertexAngle);
%                 
%                 % //////////////////// DEBUG PLOTS ////////////////////////
%                 if vertexData(i).isOnRight > 0
%                     scatter3(gca,vertex_planar(1),vertex_planar(2),vertex_planar(3),'ro');
%                 else
%                     scatter3(gca,vertex_planar(1),vertex_planar(2),vertex_planar(3),'go');
%                 end
%             end
%             % All the vertex data is now held in the structure 'vertexData'
%             candidateAngles = [vertexData(:).angle]
%             [~,rightInd]  = min(candidateAngles);
%             rightVertex = vertexData(rightInd);
%             [~,leftInd] = max(candidateAngles);
%             leftVertex = vertexData(leftInd);
%             
%             
%             % CALCULATE THE VO PARAMETERS
%             leftEdgeAngle  = leftVertex.angle   %- %asin(r_a/leftVertex.distance);
%             rightEdgeAngle = rightVertex.angle  %+ %asin(r_a/rightVertex.distance);
%             
%             % AUGMENT BY THE MINIMUM SEPARATION RADIUS
%             %             leftVertex.point + leftVertex.perp
%             
%             
%             equivalentOpenAngle  = rightEdgeAngle - leftEdgeAngle;
%             unit_leadingTangent  = OMAS_geometry.rodriguesRotation(unit_reference,planeNormal,rightEdgeAngle);
%             unit_trailingTangent = OMAS_geometry.rodriguesRotation(unit_reference,planeNormal,leftEdgeAngle);
%             
%             % /////////////// VO IS DEFINED, COMPUTE LOGIC ////////////////
%             
%             % DETERMINE WHERE Va LIES
%             if dot(unit_reference,unit_leadingTangent) > dot(unit_reference,unit_trailingTangent)
%                 isVaLeading = 1;
%             else
%                 isVaLeading = 0;
%             end
%             % DETERMINE IF Va IS INSIDE THIS VO
%             isVaInsideCone = 0;
%             VOtolerance = 1E-8;
%             
%             candidateVector = unit_reference - VOorigin;
%             
%             % EVALUATE THE PROJECTION AGAINST LEADING EDGE FIRST
%             VOprojection   = norm(candidateVector)*cos(abs(leftEdgeAngle));
%             candProjection = unit_lambda_planar'*candidateVector;
%             projDiff = (candProjection - VOprojection);
%             if projDiff > VOtolerance && isVaLeading
%                 isVaInsideCone = 1;
%             end
%             % EVALUATE THE PROJECTION AGAINST TRAILING EDGE
%             VOprojection   = norm(candidateVector)*cos(abs(rightEdgeAngle));
%             candProjection = unit_lambda_planar'*candidateVector;
%             projDiff = (candProjection - VOprojection);
%             if projDiff > VOtolerance && ~isVaLeading
%                 isVaInsideCone = 1;
%             end
%             
%             effectiveRadii = 5;
%             
            % DEFINE THE VO STRUCTURE
            VO.apex                     = origin; % VOorigin;
            VO.axisUnit                 = unit_lambda_planar;
            VO.axisLength               = norm_lambda_planar;
            VO.openAngle                = equivalentOpenAngle;
            VO.leadingEdgeUnit          = unit_leadingTangent;
            VO.trailingEdgeUnit         = unit_trailingTangent;
            VO.isVaLeading              = isVaLeading;
            VO.isVaInsideCone           = isVaInsideCone;
            VO.truncationTau            = tau;
            VO.truncationCircleCenter   = (lambda_planar)/tau + VOorigin;
            VO.truncationCircleRadius   = (effectiveRadii + r_a)/tau;
%             
%             % DEBUG PLOT
%             CCorigin = zeros(3,1);
%             scale = 10;
%             q.AutoScaleFactor = 1;
%             q = quiver3(gca,CCorigin(1),CCorigin(2),CCorigin(3),...
%                 VO.leadingEdgeUnit(1)*scale,VO.leadingEdgeUnit(2)*scale,VO.leadingEdgeUnit(3)*scale,'g');
%             q.AutoScaleFactor = 1;
%             q = quiver3(gca,CCorigin(1),CCorigin(2),CCorigin(3),...
%                 VO.trailingEdgeUnit(1)*scale,VO.trailingEdgeUnit(2)*scale,VO.trailingEdgeUnit(3)*scale,'r');
%             q.AutoScaleFactor = 1;
%             
%             % CONVERT TO 2D
%             VO.apex = VO.apex(1:2,1);
%             VO.axisUnit = OMAS_geometry.unit(VO.axisUnit(1:2,1));
%             VO.leadingEdgeUnit  = OMAS_geometry.unit(VO.leadingEdgeUnit(1:2));
%             VO.trailingEdgeUnit = OMAS_geometry.unit(VO.trailingEdgeUnit(1:2));
%             VO.truncationCircleCenter = VO.truncationCircleCenter(1:2);
        end    
        % ASSEMBLE THE 2D VELOCITY OBSTACLE (VO)
        function [VO] = define2DVelocityObstacle(obj,p_a,v_a,r_a,p_b,v_b,r_b,tau,plotOn)
            % This function assembles the standard velocity obstacle
            % in 2D. 
            
            % MAP THE 2D INPUTS TO 3D 
            p_a = [p_a;0]; 
            v_a = [v_a;0];
            p_b = [p_b;0]; 
            v_b = [v_b;0];
            % CALL THE 3D VO GENERATION FUNCTION
            VO = obj.define3DVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau,plotOn);
            % ALTER THE VO PARAMETERS 
            VO.apex = VO.apex(1:2,1);
            VO.axisUnit = VO.axisUnit(1:2,1);
            VO.leadingEdgeUnit = VO.leadingEdgeUnit(1:2,1);
            VO.trailingEdgeUnit = VO.trailingEdgeUnit(1:2,1);
            VO.truncationCircleCenter = VO.truncationCircleCenter(1:2,1);
        end   
         % CLEAR PATH STRATEGY
        function [optimalVelocity] = strategy_clearPath(obj,v_a,desiredVelocity,VOset,plotOn)
            % This function computes the optimal avoidance velocity using
            % the 'clear path' method of calculating the closest point to
            % the desired velocity on the surface of the VO set.
            
            % INPUT HANDLING
            if numel(VOset) == 0
               optimalVelocity = desiredVelocity;
               return
            end
            
            p_a = obj.localState(1:2,1);                                   % For plotting
            
            % ////////////// BUILD THE PROJECTION SET /////////////////////
            % We build a list of projection points, of closest proximity to
            % the desired velocity. There will be two projection points per
            % VO. 
            projectionPoints = zeros(2,2*numel(VOset));
            isOnRayPoints    = ones(1,2*numel(VOset));
            a = 0;  
            for VOnumA = 1:numel(VOset)
                % THE FIRST VERTEX EDGE
                [projections(:,1),isOnRay(1)] = obj.pointProjectionToRay(desiredVelocity,VOset(VOnumA).apex,VOset(VOnumA).leadingEdgeUnit);
                % THE SECOND VERTEX EDGE
                [projections(:,2),isOnRay(2)] = obj.pointProjectionToRay(desiredVelocity,VOset(VOnumA).apex,VOset(VOnumA).trailingEdgeUnit);
                % COLLECT THE PROJECTIONS POINTS
                % The projections of 'v_a' on both the leadingEdgeUnit, trailingEdgeUnit
                isOnRayPoints((1 + a*VOnumA):(2 + a*VOnumA)) = isOnRay;          % CONCATINATE THE IS ON RAY 
                projectionPoints(:,(1 + a*VOnumA):(2 + a*VOnumA)) = projections; % STORE ALL PROJECTION POINTS
                a = a + 1;
            end

            % /////////// BUILD THE INTERSECTION POINT SET ////////////////
            % GET THE INTERSECTIONS BETWEEN TWO SETS OF LEADING & TRAILING
            % EDGES
            VOsum = numel(VOset);
            intersectionFlags  = ones(1,4*VOsum*(VOsum-1)/2);
            intersectionPoints = zeros(2,4*VOsum*(VOsum-1)/2);
            a = 0;
            for VOnum_i = 1:numel(VOset)
                for VOnum_j = 1:numel(VOset)
                    if VOnum_i == VOnum_j 
                       continue % Skip self comparison (also omits singular VO condition) 
                    end
                    pIntersect = zeros(2,4);
                    % LEADING - LEADING
                    [pIntersect(:,1),validIntersect(1)] = obj.twoRayIntersection2D(...
                        VOset(VOnum_i).apex,...
                        VOset(VOnum_i).leadingEdgeUnit,...
                        VOset(VOnum_j).apex,...
                        VOset(VOnum_j).leadingEdgeUnit);
                    % LEADING - TRAILING
                    [pIntersect(:,2),validIntersect(2)] = obj.twoRayIntersection2D(...
                        VOset(VOnum_i).apex,...
                        VOset(VOnum_i).leadingEdgeUnit,...
                        VOset(VOnum_j).apex,...
                        VOset(VOnum_j).trailingEdgeUnit);                    
                    % TRAILING - LEADING
                    [pIntersect(:,3),validIntersect(3)] = obj.twoRayIntersection2D(...
                        VOset(VOnum_i).apex,...
                        VOset(VOnum_i).trailingEdgeUnit,...
                        VOset(VOnum_j).apex,...
                        VOset(VOnum_j).leadingEdgeUnit);
                    % TRAILING - TRAILING
                    [pIntersect(:,4),validIntersect(4)] = obj.twoRayIntersection2D(...
                        VOset(VOnum_i).apex,...
                        VOset(VOnum_i).trailingEdgeUnit,...
                        VOset(VOnum_j).apex,...
                        VOset(VOnum_j).trailingEdgeUnit);        
                    
                    % There are four intersections per pair of VO. 

                    % RETAIN THE POINTS & FLAGS
                    intersectionFlags(:,(1 + 4*a):4*(1 + a))  = validIntersect; % If the corresponding point was a valid intersection
                    intersectionPoints(:,(1 + 4*a):4*(1 + a)) = pIntersect;     % The intersection point array
                    a = a + 1;
                end
            end
            
            % ASSSESS THE COLLECTIVE POINT SET AGAINSTS THE VO SET
            % All valid projections and intersection must be compared
            % against thw VO set.
            
            % OMIT NON-VALID PROJECTIONS
            validProjectionPoints = projectionPoints(:,(isOnRayPoints == 1));    % Get only the projections where the points are on rays            
            % REMOVE ANY NON-INTERSECTIONS
            validIntersectionPoints = intersectionPoints(:,(intersectionFlags == 1)); % Are valid intersections

            % CONSIDER THE CURRENT VELOCITY IN THE CANDIDATE SET
            collectivePoints = [desiredVelocity,validProjectionPoints,validIntersectionPoints]; % <<< TO BE CONFIRMED
            collectivePoints = unique(collectivePoints','rows');           % Remove repeat candidates
            collectivePoints = collectivePoints';
            
            % ///////// CHECK EACH POINT AGAINST THE VO SET ///////////////
            VOflagVector = zeros(1,size(collectivePoints,2));
            for candidate = 1:size(collectivePoints,2)
                for VOnum_i = 1:numel(VOset)
                    isInsideVOFlag = obj.isInsideVO(VOset(VOnum_i).apex,...
                                                    VOset(VOnum_i).axisUnit,...
                                                    VOset(VOnum_i).leadingEdgeUnit,...
                                                    VOset(VOnum_i).trailingEdgeUnit,...
                                                    collectivePoints(:,candidate));
                    % DETERMINE WHETHER THE POINT BELONGS TO ANY VO
%                     if VOflagVector(candidate) || obj.isInsideVO(collectivePoints(:,candidate),VOset(VOnum_i))
                    if VOflagVector(candidate) || isInsideVOFlag                   
                    	VOflagVector(candidate) = 1;
                    end
                end
            end
            
            % REMOVE THE VO-INVALIDATED CANDIDATE POINTS
            candidatesOutsideVO = collectivePoints(:,VOflagVector ~= 1);
            
            % ///// CHOOSE OPTIMAL VELOCITY FROM THE CANDIDATE POINTS /////
            compareVelocity = desiredVelocity;
            optimalVelocity = zeros(2,1);
            
            % //////////////// THE MINIMUM DIFFERENCE /////////////////////
            optimalMetricDistance = inf;  % Metric of optimality
            % ABSOLUTE MAGNITUDE COMPARISON
            if size(candidatesOutsideVO,2) > 0 
                % ASSESS VELOCITIES AGAINST THE DESIRED VELOCITY
                for k = 1:size(candidatesOutsideVO,2)
                    dis = norm(candidatesOutsideVO(:,k) - compareVelocity);
                    if dis < optimalMetricDistance 
                        optimalVelocity = candidatesOutsideVO(:,k);
                        optimalMetricDistance = dis;
                    end
                end
            elseif isempty(candidatesOutsideVO)
                % IN THE EVENT THERE ARE NO VALID VELOCITIES
                warning('There is no feasible velocity!');
                optimalVelocity = zeros(2,1);                
            end

            % /////////// THE SMALLEST HEADING DIFFERENCE /////////////////
%             unit_compareVelocity = OMAS_geometry.unit(compareVelocity);
%             optimalHeadingChange = 0;
%             % HEADING COMPARISON
%             if size(candidatesOutsideVO,2) > 0 
%                 % ASSESS VELOCITIES AGAINST THE DESIRED VELOCITY
%                 for k = 1:size(candidatesOutsideVO,2)
%                     % MEASURE THE ANGLE BETWEEN CANDIDATES
%                     dotProduct = abs(dot(unit_compareVelocity,candidatesOutsideVO(:,k)));
%                     
%                     if dotProduct > optimalHeadingChange
%                         optimalVelocity = candidatesOutsideVO(:,k);
%                         optimalHeadingChange = dotProduct;                 % Set as new target to beat
%                     end
%                 end
%             elseif isempty(candidatesOutsideVO)
%                 % IN THE EVENT THERE ARE NO VALID VELOCITIES
%                 warning('There is no feasible velocity!');
%                 optimalVelocity = zeros(2,1);                
%             end 
            
            % //////////////////  PLOT THE SCENE //////////////////////////
            if obj.objectID == 1 && plotOn
                hold on; grid on;
                origin = zeros(3,1);    % DRAW RELATIVE TO ORIGIN
                
                for proj_i = 1:size(projectionPoints,2)
                    plotProj = projectionPoints(:,proj_i);
                    scatter(gca,plotProj(1),plotProj(2),'r');
                end
                
                for inter_j = 1:size(intersectionPoints,2)
                    plotProj = intersectionPoints(:,inter_j);
                    scatter(gca,plotProj(1),plotProj(2),'b');
                end
                
                % FILL IN VALID POINTS
                for validPoint = 1:size(candidatesOutsideVO,2)
                   scatter(gca,candidatesOutsideVO(1,validPoint),candidatesOutsideVO(2,validPoint),'b','filled'); 
                end
                
                % OUR VELOCITY
                q = quiver(gca,origin(1),origin(2),v_a(1),v_a(2),'r');
                q.AutoScaleFactor = 1;
                
                % DESIRED VELOCITY
                
                q = quiver(gca,origin(1),origin(2),desiredVelocity(1),desiredVelocity(2),'b');
                q.AutoScaleFactor = 1;
                scatter(gca,desiredVelocity(1),desiredVelocity(2),'b');
                
                % THE OPTIMAL VELOCITY
                q = quiver(gca,origin(1),origin(2),optimalVelocity(1),optimalVelocity(2),'m');
                q.AutoScaleFactor = 1;
                scatter(gca,optimalVelocity(1),optimalVelocity(2),'m');
                
                for VOnum_i = 1:numel(VOset)
                    leadingEdgeVector  = VOset(VOnum_i).axisLength*VOset(VOnum_i).leadingEdgeUnit;
                    trailingEdgeVector = VOset(VOnum_i).axisLength*VOset(VOnum_i).trailingEdgeUnit;
                    apexPosition = VOset(VOnum_i).apex;
                    q = quiver(gca,apexPosition(1),apexPosition(2),leadingEdgeVector(1),leadingEdgeVector(2),'b');
                    q.AutoScaleFactor = 1;
                    q = quiver(gca,apexPosition(1),apexPosition(2),trailingEdgeVector(1),trailingEdgeVector(2),'r');
                    q.AutoScaleFactor = 1;
                end
            end
        end       
    end   
    % //////////////////////////// UTILITIES //////////////////////////////
    methods (Static) 
        % TWO LINE INTERSECTIONS (FROM agent_VO)
        function [p_inter,isSuccessful] = twoRayIntersection2D(P1,dP1,P2,dP2)
            % Find the intersection point between two 2D vectors. This
            % function isnt' interested if vertices are infront or behind
            % of the starting point.
            % INPUTS:
            % - P1,P2   - The ray defining points.
            % - dP1,dP2 - The ray unit directions.
            % OUTPUTS:
            % - p_inter - The 2D intersection point.
            
            assert(numel(P1) == 2,'Input must be 2D');
            assert(numel(P2) == 2,'Input must be 2D');
            
            % SOME SUFFICIENTLY SMALL VALUE FOR INTERSECTION
            isSuccessful = logical(false);   % Default to no intersection
            p_inter = NaN(2,1);              % Default to no intersection
            
            % THE 2D DETERMININANT
            div = dP1(2)*dP2(1) - dP1(1)*dP2(2);
            if div == 0
                disp('Lines are parallel');
                return % Lines are parallel
            end
            
            % THE SCALAR PROJECTIONS
            mua = (dP2(1)*(P2(2) - P1(2)) + dP2(2)*(P1(1) - P2(1))) / div;
            mub = (dP1(1)*(P2(2) - P1(2)) + dP1(2)*(P1(1) - P2(1))) / div;
            
            % POINTS MUST BE THE RESULT OF A POSITIVE INCREMENT OF THE VECTOR GRADIENT
            % (i.e, in the correct direction)
            if mua < 0 || mub < 0   % Intersections only occur in the direction of the vector
                return              % Lines do not intersect
            end
            
            % THE INTERSECTION POINT
            p_inter = P1 + mua*dP1;
            isSuccessful = logical(true);
        end
        % GET THE POINT PROJECTION ON RAY
        function [projectedPoint,isOnTheRay] = pointProjectionToRay(p,p0,v0)
            % INPUTS:
            % p  - Is the point to be projected.
            % p0,v0 - The line defining points
            % OUTPUTS:
            %
            projectedPoint = v0*v0'/(v0'*v0)*(p - p0) + p0;
            
            if v0'*(projectedPoint - p0)>0 % if on the ray
                isOnTheRay = logical(true);
            else
                isOnTheRay = logical(false);
            end
        end
        % CHECK THE AGENT VELOCITY IS INSIDE AN (ASYMMETRICAL) VO
        function [flag] = isInsideVO(VOapex,VOaxisUnit,VOedgeA,VOedgeB,point)
            % determine if the point p is inside the given VO
            % angle<halfOpenAngle; here numericalTolerance is an error tolarance
                        
            flag = 0;
            VOtolerance = 1E-8;
            
            % Vector to point
            candidateVector = point - VOapex;
            unit_candidate = OMAS_geometry.unit(candidateVector);
            
            % Candidate projection
            normVector = cross([VOaxisUnit;0],[VOedgeA;0]);
            candNorm = cross([VOaxisUnit;0],[unit_candidate;0]);
            
            % Direction (+ve) if comparing unit_candidate and edge A are on
            % the same side of the VO axis
            direction = sign(dot(normVector,candNorm));
            
            candVerProj = dot(VOaxisUnit,unit_candidate);
            
            % Determine which edge the candidate vector  should be compared to
            if direction > 0
            	VOprojection = dot(VOaxisUnit,VOedgeA);
            else
                VOprojection = dot(VOaxisUnit,VOedgeB);
            end
            
            projDiff = candVerProj - VOprojection;
            if projDiff > VOtolerance
                flag = 1;
            end
            
%             f = figure(2);
%             hold on;
%             q = quiver(gca,0,0,VOaxisUnit(1),VOaxisUnit(2),'k');
%             q.AutoScaleFactor = 1;
% 
%             q = quiver(gca,0,0,VOedgeA(1),VOedgeA(2),'r');
%             q.AutoScaleFactor = 1;
%             q = quiver(gca,0,0,VOedgeB(1),VOedgeB(2),'g');
%             q.AutoScaleFactor = 1;
%             
%             q = quiver(gca,0,0,unit_candidate(1),unit_candidate(2),'b');
%             q.AutoScaleFactor = 1;    
%             hold off;
            
%             % Evaluate based on symmetry of the VO
%             if isSymmetric
%                 VOprojection = norm(candidateVector)*cos(VOopenAngle/2);
%                 candProjection = VOaxisUnit'*candidateVector;
%                 projDiff = (candProjection - VOprojection);
%                 if projDiff > VOtolerance   
%                     flag = 1;
%                 end
%             else
% %                 VOprojection = norm(candidateVector)*cos(VO.openAngle/2);
%             end

        end
        % GET THE UNIT VELOCITY WITH NAN PROTECTION
        function [unitVelocity,mag] = nullVelocityCheck(velocity)
            % This function is designed to prevent the NaN propagation
            % occurring when the a zero magnititude velocity is used to
            % represent the direction of the obstacle.
            
            tolerance = 1E-8;
            
            % NULL VELOCITY CATCHA
            if sum(velocity) < tolerance
                unitVelocity = zeros(size(velocity));
                mag = 0;
            else
                mag = norm(velocity);
                unitVelocity = velocity/mag;
            end
        end
    end
    % //////////////////// CONTROLLER/ PLANT OPERATIONS ///////////////////
    methods
        % CONTROLLER
        function [obj] = controller(obj,dt,desiredVelocity)
            % This function computes the agents change in state as a result
            % of the desired velocity vector
            
            % SANITY CHECK OF THE BOUNDED VELOCITY
            [unitDirection,desiredSpeed] = obj.nullVelocityCheck(desiredVelocity);   
            
            % APPLY SPEED CONSTRAINT
            if desiredSpeed > obj.maxSpeed
                desiredSpeed = obj.maxSpeed; 
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dHeading] = obj.GetVectorHeadingAngles([1;0;0],[unitDirection;0]); % Relative heading angles   
            omega = -dHeading/dt;
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:3),[desiredSpeed;0],omega);
            obj.localState(1:3) = obj.localState(1:3) + dt*dX;
            obj.localState(4:6) = dX;

            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES /////////////// 
            if obj.VIRTUAL.idleStatus
                obj.localState(4:5) = zeros(2,1);
            end
            % THE 2D VELOCITY 
            [obj] = obj.updateGlobalProperties_2DVelocities(dt,obj.localState);
        end
    end
end
% AGENT STATE VECTOR [x;y;psi;xdot;ydot;psidot]