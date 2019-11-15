%% THE RECIPROCAL VELOCITY OBSTACLE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This agent class provides support for the RVO2 library for ORCA-based
% Reciprocal Velocity Obstacle collision avoidance.

% Author: James A. Douthwaite

classdef agent_2D_RVO2 < agent_2D_RVO
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % AGENT PARAMETERS
%         radius = 2;          % Declare (independant of the VO varients)
%         neighbourDist = 15;  %(s)
%         maxNeighbours = 10;  
        timeHorizon     = 10;
        timeHorizonObst = 10;
        RVO_EPSILON     = 1E-5;
        
        % PARAMETERS UNIQUE TO RVO2
%         timeHorizon     = 10;
%         timeHorizonObst = 10;
                
        % ORCA-LINE CONTAINERS
        orcaLines = struct('point',[],'direction',[]);
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_2D_RVO2(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_RVO(varargin); 

            % SENSORS
            % Sensor properties will be inherited from 'agent_2D_VO' and 
            % 'agent_VO' unless explictly overridden here.

            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin); 
            % /////////////////////////////////////////////////////////////
        end
    end
    % //////////////////////// OPENMAS FUNCTIONS ///////////////////////////
    methods 
        % MAPPING OF OPENMAS TO THE RVO LIBRARY
        function [headingVector,speed] = GetAvoidanceCorrection(obj,dt,desiredVelocity,visualiseProblem)
            % This function computes the mapping of the simulation inputs
            % for the RVO library functions and returns the calculated
            % optimal velocity.
            
            % AGENT KNOWLEDGE (2D)
            [p_i,v_i,r_i] = obj.GetAgentMeasurements();
            
            % PLOT INPUT VARIABLES
            inputDebugPlot = 0;
            if inputDebugPlot 
                figureHandle = figure(1);
                hold on; grid on; box on;
                axis equal;
                [figureHandle] = obj.getObjectScene(knownObstacles,figureHandle);           % Plot the objects
            end
            
            % GET AGENT DATA
            agentIDs = [obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.agent).objectID];
            
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            agentNeighbours = cell(numel(agentIDs),1); 
            agentNeighbourLength = 0;
            for item = 1:numel(agentIDs)
                % Get object data from memory structure
                p_j = obj.GetLastMeasurementByID(agentIDs(item),'position');
                v_j = obj.GetLastMeasurementByID(agentIDs(item),'velocity');
                r_j = obj.GetLastMeasurementByID(agentIDs(item),'radius');
                
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < obj.neighbourDist;       % [CONFIRMED] 
                neighbourConditionC = ~any(isnan(v_j));                    % Wait for a valid velocity reading
                if ~neighbourConditionB || ~neighbourConditionC
                    continue
                end
            
                p_j = p_j + p_i; 
                v_j = v_j + v_i;  
                
                % BUILD THE RVO AGENT DESCRIPTION  
                agentRef = struct('position',p_j,...
                                  'velocity',v_j,...
                                  'unitDir',obj.nullVelocityCheck(v_j),...
                                  'radius',r_j,...
                                  'id',agentIDs(item));
                              
                % BUILD LIST OF AGENTS
%                 agentNeighbours = horzcat(agentNeighbours,agentRef);
                agentNeighbours{item} = agentRef;
                agentNeighbourLength  = agentNeighbourLength + 1;
            end
            % REMOVE EMPTY CELLS
            agentNeighbours = agentNeighbours(~cellfun('isempty',agentNeighbours));  
            
            % OBSTACLE DATA
            obstacleIDs  = [obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.obstacle).objectID];
            
            % BUILD THE REQUIRED OBSTACLE STRUCTURES
            obstacleNeighbours = cell(numel(obstacleIDs),1);
            for item = 1:numel(obstacleIDs)
                % Get object data from memory structure
                p_j = obj.GetLastMeasurementByID(obstacleIDs(item),'position');
                v_j = obj.GetLastMeasurementByID(obstacleIDs(item),'velocity');
              
                % PROCESS THE OBSTACLE'S VERTICES
                obstacleRef.position = p_j + p_i;
                obstacleRef.velocity = v_j + v_i;
%               obstacleRef.geometry.vertices = obstacleRef.geometry.vertices + [agentPosition;0]';
                g_j = obj.GetLastMeasurementByID(obstacleIDs(item),'geometry');
                
                % [TO-DO] PASS THE OBSTACLE GEOMETRY 
                obstacleNeighbours{item} = obj.createRVOobstacleVertexList(g_j);
            end
            % REMOVE EMPTY CELLS
            obstacleNeighbours = obstacleNeighbours(~cellfun('isempty',obstacleNeighbours));
                       
            % COMPUTE THE AGENT'S NEW VELOCITY
            % This calls the libraries 'computeNewVelocity' function
            [obj,avoidanceVelocity] = obj.computeNewVelocity(dt,desiredVelocity,agentNeighbours,obstacleNeighbours);
        
            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0];  % Retain previous heading
            end
        
        end
    end
    %% FUNCTIONS IN THE "agent.cpp" FILE
    methods
        % COMPUTE THE AGENT VELOCITY
        function [obj,avoidanceVelocity] = computeNewVelocity(obj,dt,prefVelocity,agentNeighbours,obstacleNeighbours)
            % This function computes the new velocity of this agent from the
            % RVO2 library.
            
            % DEFINE THE AGENT PARAMETERS
            [p_a,v_a,r_a] = obj.GetAgentMeasurements();

            % MAP THE PREFERRED VELOCITY TO NEW FRAME
            prefVelocity = -prefVelocity;
                     
            % RESET THE OCRA-LINES VECTOR FROM THE LAST TIMESTEP
            obj.orcaLines = [];
              
            % //////////// GET THE OBSTACLE ORCA LINES ////////////////////
            [obst_orcaLines] = obj.getObstacleORCAlines(p_a,v_a,r_a,obj.timeHorizonObst,obstacleNeighbours);
            obj.orcaLines = [obj.orcaLines,obst_orcaLines];
            numObstLines = numel(obj.orcaLines);
            
            % ////////////// GET THE AGENT ORCA LINES /////////////////////
            [agnt_orcaLines] = obj.getAgentORCAlines(dt,p_a,v_a,r_a,obj.timeHorizon,agentNeighbours);
            obj.orcaLines = [obj.orcaLines;agnt_orcaLines];           
            
            % //////////////// DEBUG PLOTS ////////////////////////////////
%             if ~isempty(obj.orcaLines)
%                 % PLOT THE ORCA-LINES
%                 h = figure(1);
%                 ax = axes(h);
%                 hold on; grid on; box on;
%                 axis equal;
% 
%                 g_a = obj.GEOMETRY.vertices + [p_a;0]';
%                 patch(ax,...
%                     'Faces',obj.GEOMETRY.faces,...
%                     'Vertices',g_a);
%                 
%                 % The objects vertices
%                 objectVertices = obstacleNeighbours{1};
%                 for i = 1:numel(objectVertices)
%                     scatter(ax,objectVertices(i).point(1),objectVertices(i).point(2),'r*'); 
%                 end  
% 
%                 for i = 1:numel(obj.orcaLines)
%                     lineVar = obj.orcaLines(i);
%                     len = 5;
%                     X = [lineVar.point(1);len*lineVar.direction(1)];
%                     Y = [lineVar.point(2);len*lineVar.direction(2)];
%                     plot(ax,X,Y,'b','marker','o','lineWidth',2);
%                 end
%                 close(h);
%             end
            
            % CALL THE SECOND LINEAR PROGRAM
            [obj,lineObjFail,optimalVelocity] = obj.linearProgram2(obj.orcaLines,obj.v_max,prefVelocity,logical(false));
            
            % IF SOME ORCA-LINES FAIL
            % This function will generate direction based ORCA lines where
            % the previous lines failed.
            if lineObjFail < numel(obj.orcaLines)
                [obj,optimalVelocity] = obj.linearProgram3(obj.orcaLines,numObstLines,lineObjFail,obj.v_max,optimalVelocity);
            end
            % REASSIGN OUTPUT VELOCITY
%             avoidanceVelocity = optimalVelocity;
            avoidanceVelocity = -optimalVelocity;
        end
        % PROCESS THE VERTEX SET
        function [vertexSet] = createRVOobstacleVertexList(obj,observedObject)
            % The obstacle consideration expects the obstacles to be of
            % polygonal form. Each obstacle is defined by a subset of
            % points that resemble its vertices.
            
            % We assume that each obstacle defines a set of obstacle
            % vertices. In the RVO2 library, the assumption is made
            % that the vertices are given in counter-clockwise order.
            
            %/*
            %* Add (polygonal) obstacles, specifying their vertices in counterclockwise
            %* order.
            %*/
            
            % THE OBSTACLES VERTICES
            observedVertices = observedObject.geometry.vertices;
            
            % GET THE VERTICES PROJECTED ONTO THE X/Y PLANE
            [planarVertices] = OMAS_graphics.geometryPlanarProjection([0;0;1],...
                                                                      observedObject.position,...
                                                                      observedVertices);
            % SORT THE VERTICES CCW ABOUT THE +ve Z-AXIS 
            [verticesCCW] = OMAS_graphics.sortVerticesCCW([0;0;1],...
                                                          observedObject.position,...
                                                          planarVertices);
            % GETTING THE 2D VERTEX SET
            if size(verticesCCW,2) > 2
               verticesCCW = verticesCCW(:,1:2);             
            end                                                    
            % BUILD THE RVO POLY-GONAL OSTACLE DESCRIPTION  
            [vertexSet] = obj.defineRVOVertexObstacleSet(verticesCCW);
        end
        % GET THE OBSTACLE ORCA-LINE VECTOR % <<<<<<<<< THE FINAL PROBLEM
        function [obst_orcalines] = getObstacleORCAlines(obj,pa,va,ra,timeHorizonObst,obstacleNeighbours)
            % This function calculates the the orca-lines for static
            % obstacles using the prinicples of the VO.
            % INPUTS:
            % obj - Agent self-reference
            % dt  - The simulation timestep
            
            % The obstacle neighbour list is composed of the discrete 
            % vertex sets. Each neighbour is a set of vertices with their
            % links to each other specified.
            
            invTimeHorizonObst = 1/timeHorizonObst;
            
            % ORCA-LINE CONTAINER
            obst_orcalines = [];
            
            if isempty(obstacleNeighbours)
                return
            end
            
            % FOR EACH NEIGHBOURING OBSTACLE
            for n = 1:numel(obstacleNeighbours)
                % THE SUBSET OF VERTICES
                obstacleVertexSet = obstacleNeighbours{n};
                % [JD] - Small notation change. The obstacleVertexSet is
                % now a list of structures that defines the relationship
                % between each vertex. 
                
                % FOR EACH VERTEX IN THE NEIGHBOUR SET
                for i = 1:numel(obstacleVertexSet)
                    % GET THE OBSTACLE
                    obstacle1 = obstacleVertexSet(i);
                    obstacle2 = obstacleVertexSet(obstacle1.nextObstacle);
                    % RELATIVE POSITION OF THE POINTS
                    relativePosition1 = obstacle1.point - pa;
                    relativePosition2 = obstacle2.point - pa;
                    
                    % Check if velocity obstacle of obstacle is already taken care of by
                    % previously constructed obstacle ORCA lines.
                    alreadyCovered = logical(false);

                    if ~isempty(obst_orcalines)
                        for j = 1:numel(obj.orcaLines)
                            conditionA = obj.det(invTimeHorizonObst*relativePosition1 - obst_orcalines(j).point, obst_orcalines(j).direction) - invTimeHorizonObst*ra >= -obj.RVO_EPSILON;
                            conditionB = obj.det(invTimeHorizonObst*relativePosition2 - obst_orcalines(j).point, obst_orcalines(j).direction) - invTimeHorizonObst*ra >= -obj.RVO_EPSILON;
                            if conditionA && conditionB
                                alreadyCovered = logical(true);
                                break;
                            end
                        end
                    end
                    % If covered, move to the next obstacle
                    if alreadyCovered
                        continue
                    end

                    % Obstacle is not yet covered, check for collisions
                    distSq1 = obj.absSq(relativePosition1);
                    distSq2 = obj.absSq(relativePosition2);

                    radiusSq = ra^2;

                    obstacleVector = obstacle2.point - obstacle1.point;
                    s = obj.dot(-relativePosition1,obstacleVector)/obj.absSq(obstacleVector);
                    distSqlineObj = obj.absSq(-relativePosition1 - s*obstacleVector);

                    % ////////////// CREATE ORCA (lineObj) ////////////////
                    lineObj = struct('point',[],'direction',[]);

                    if s < 0 && distSq1 <= radiusSq
                        % Collision with left vertex. Ignore if non-convex
                        if obstacle1.isConvex
                            lineObj.point = [0;0];
                            direction = [-relativePosition1(2);relativePosition1(1)];  % [y;x]
                            lineObj.direction = obj.normalise(direction);
                            obst_orcalines = [obst_orcalines;lineObj];
                        end
                        continue
                    elseif s > 1 && distSq2 <= radiusSq
                        % Collision with right vertex. Ignore if non-convex or
                        % if it will be taken care of by neighbouring obstacle.
                        if obstacle2.isConvex && obj.det(relativePosition2,obstacle2.unitDir) >= 0
                            lineObj.point = [0;0];
                            direction = [-relativePosition2(2);relativePosition2(1)];
                            lineObj.direction = obj.normalise(direction);
                            obst_orcalines = [obst_orcalines;lineObj];
                        end
                        continue
                    elseif s >= 0 && s < 1 && distSqlineObj <= radiusSq
                        % Collision with obstacle segment
                        lineObj.point = [0;0];
                        lineObj.direction = -obstacle1.unitDir;
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue
                    end

%                     /*
%                      * No collision.
%                      * Compute legs. When obliquely viewed, both legs can come from a single
%                      * vertex. Legs extend cut-off line when nonconvex vertex.
%                      */

                    if s < 0 && distSqlineObj <= radiusSq
                        % Obstacle viewed obliquely so that left vertex defines
                        % velocity obstacle
                        if ~obstacle1.isConvex
                            continue;   % Ignore obstacle
                        end
                        obstacle2 = obstacle1;

                        leg1 = sqrt(distSq1 - radiusSq);
                        leftLegDirection  = [relativePosition1(1)*leg1 - relativePosition1(2)*ra;  relativePosition1(1)*ra + relativePosition1(2)*leg1] / distSq1;
                        rightLegDirection = [relativePosition1(1)*leg1 + relativePosition1(2)*ra; -relativePosition1(1)*ra + relativePosition1(2)*leg1] / distSq1; % Assuming [x;y]
                    elseif s > 1 && distSqlineObj <= radiusSq
                        % Obstacle viewed obliquely so that right vertex defines
                        % velocity obstacle.
                        if ~obstacle2.isConvex
                            % Ignore obstacle
                            continue
                        end
                        obstacle1 = obstacle2;

                        leg2 = sqrt(distSq2 - radiusSq);
                        leftLegDirection  = [relativePosition2(1)*leg2 - relativePosition2(2)*ra;  relativePosition2(1)*ra + relativePosition2(2)*leg2] / distSq2;
                        rightLegDirection = [relativePosition2(1)*leg2 + relativePosition2(2)*ra; -relativePosition2(1)*ra + relativePosition2(2)*leg2] / distSq2; % Assuming [x;y]
                    else
                        % Unusual situation
                        if obstacle1.isConvex
                            leg1 = sqrt(distSq1 - radiusSq);
                            leftLegDirection = [relativePosition1(1)*leg1 - relativePosition1(2)*ra; relativePosition1(1)*ra + relativePosition1(2)*leg1] / distSq1;
                        else
                            % Left vertex non-convex; left leg extends cut-off
                            % lineObj
                            leftLegDirection = -obstacle1.unitDir;
                        end

                        if obstacle2.isConvex
                            leg2 = sqrt(distSq2 - radiusSq);
                            rightLegDirection = [relativePosition2(1)*leg2 + relativePosition2(2)*ra; -relativePosition2(1)*ra + relativePosition2(2)*leg2] / distSq2;
                        else
                            % Right vertex non-convex; right leg extends cut-off
                            % lineObj
                            rightLegDirection = obstacle1.unitDir;
                        end
                    end

                    % Legs can never point into neighbouring edge when convex
                    % vertex, take cutoff-lineObj of neighbouring edge instead. If
                    % velocity projected on "foreign" leg, no constraint is
                    % added.

                    leftNeighbour = obstacleVertexSet(obstacle1.prevObstacle);
                    isLeftLegForeign  = logical(false);
                    isRightLegForeign = logical(false);

                    if obstacle1.isConvex && obj.det(leftLegDirection,-leftNeighbour.unitDir) >= 0
                        % Left leg points into obstacle.
                        leftLegDirection = -leftNeighbour.unitDir;
                        isLeftLegForeign = logical(true);
                    end
                    if obstacle2.isConvex && obj.det(rightLegDirection, obstacle2.unitDir) <= 0
                        % Right leg points into obstacle.
                        rightLegDirection = obstacle2.unitDir;
                        isRightLegForeign = logical(true);
                    end
                    % COMPUTE CUT-OFF CENTERS
                    leftCutoff  = invTimeHorizonObst*(obstacle1.point - pa);
                    rightCutoff = invTimeHorizonObst*(obstacle2.point - pa);
                    cutoffVec   = rightCutoff - leftCutoff;

                    % Project current velocity on velocity obstacle.

                    % Check if current velocity is projected on cutoff circles
                    if obstacle1.id == obstacle2.id
                        t = 0.5;
                    else
                        t = obj.dot((va - leftCutoff),cutoffVec) / obj.absSq(cutoffVec);
                    end

                    tLeft  = obj.dot((va - leftCutoff),leftLegDirection);
                    tRight = obj.dot((va - rightCutoff),rightLegDirection);

                    if ((t < 0 && tLeft < 0) || (obstacle1.id == obstacle2.id && tLeft < 0 && tRight < 0))
                        % Project on left cut-off circle.
                        unitW = obj.normalise(va - leftCutoff);

                        lineObj.direction = [unitW(2); -unitW(1)];
                        lineObj.point = leftCutoff + ra*invTimeHorizonObst*unitW;
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue
                    elseif (t > 1 && tRight < 0)
                        % Project on right cut-off circle.
                        unitW = obj.normalise(va - rightCutoff);

                        lineObj.direction = [unitW(2); -unitW(1)];
                        lineObj.point = rightCutoff + ra*invTimeHorizonObst*unitW;
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue;
                    end

                    % Project on left leg, right leg, or cut-off lineObj, whichever is closest
                    % to velocity.
                    % CALCULATE distSqCutoff
                    if t < 0 || t > 1 || obstacle1.id == obstacle2.id
                        distSqCutOff = inf;
                    else
                        distSqCutOff = obj.absSq(va - (leftCutoff + t*cutoffVec));
                    end

                    % CALCULATE distSqLeft
                    if (tLeft < 0)
                        distSqLeft = inf;
                    else
                        distSqLeft = obj.absSq(va - (leftCutoff + tLeft*leftLegDirection));
                    end
                    % CALCULATE distSqRight
                    if (tRight < 0)
                        distSqRight = inf;
                    else
                        distSqRight = obj.absSq(va - (rightCutoff + tRight*rightLegDirection));
                    end

                    if distSqCutOff <= distSqLeft && distSqCutOff <= distSqRight
                        % Project on cut-off lineObj.
                        lineObj.direction = -obstacle1.unitDir;
                        lineObj.point = leftCutoff + ra*invTimeHorizonObst*[-lineObj.direction(2); lineObj.direction(1)];
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue;
                    elseif distSqLeft <= distSqRight
                        % Project on left leg.
                        if isLeftLegForeign
                            continue
                        end
                        lineObj.direction = leftLegDirection;
                        lineObj.point = leftCutoff + ra*invTimeHorizonObst*[-lineObj.direction(2); lineObj.direction(1)];
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue;
                    else
                        % Project on right leg. */
                        if isRightLegForeign
                            continue;
                        end
                        lineObj.direction = -rightLegDirection;
                        lineObj.point = rightCutoff + ra*invTimeHorizonObst*[-lineObj.direction(2); lineObj.direction(1)];
                        obst_orcalines = [obst_orcalines;lineObj];
                        continue;
                    end
                    
                end
            end
        end
        % GET THE AGENT ORCA-LINE VECTOR
        function [agnt_orcalines] = getAgentORCAlines(obj,dt,pa,va,ra,timehorizon,agentNeighbours)
            % This function computes the orca lines using the principles of
            % the RVO method.
            % INPUTS:
            % obj - Agent reference, for use of local functions.
            % dt  - The simulation timestep.
            % pa  - The agents global position
            % va  - The agents global velocity
            % ra  - The agents characteristic radius
            % timeHorizon - The time for which avoidance should be assured.
            % agentNeighbours - The list of considered agent-neighbours.
            
            % ORCA-LINE CONTAINER
            agnt_orcalines = [];
            
            if isempty(agentNeighbours)
                return
            end
            
            % BEGIN CREATING ORCA LINES FOR THE AGENTS
            invTimeHorizon  = 1/timehorizon;
            
            for i = 1:numel(agentNeighbours)
                % CALCULATE THE OTHER AGENT PROPERTIES
%                 other = agentNeighbours(i);
                %                     relativePosition = other.position - pa;
                %                     relativeVelocity = va - other.velocity;
%                 agentNeighbours{i}
                relativePosition = -(agentNeighbours{i}.position - pa);
                relativeVelocity = -(va - agentNeighbours{i}.velocity);
                
                distSq = obj.dot(relativePosition,relativePosition);
                combinedRadius = ra + agentNeighbours{i}.radius;
                combinedRadiusSq = combinedRadius^2;
                
                % CREATE ORCA lineObj
                lineObj = struct('point',[],'direction',[]);
                
                % NO COLLISION
                if distSq > combinedRadiusSq
                    w = relativeVelocity - invTimeHorizon*relativePosition; % Vector from cutoff center to relative velocity
                    wLengthSq   = obj.dot(w,w); % Square length
                    dotProduct1 = obj.dot(w,relativePosition);
                    
                    if dotProduct1 < 0 && dotProduct1^2 > combinedRadiusSq*wLengthSq
                        % Project on cut-off circle
                        wLength = sqrt(wLengthSq);
                        unitW = w/wLength;
                        
                        lineObj.direction = [unitW(2); -unitW(1)];
                        u = (combinedRadius*invTimeHorizon - wLength)*unitW;
                        
                    else
                        % Project on legs
                        leg = sqrt(distSq - combinedRadiusSq);
                        
                        if obj.det(relativePosition, w) > 0
                            % Project on left leg
                            lineObj.direction = ([relativePosition(1)*leg - relativePosition(2)*combinedRadius;
                                relativePosition(1)*combinedRadius + relativePosition(2)*leg]) / distSq;
                        else
                            % Project on right leg
                            lineObj.direction = -([relativePosition(1)*leg + relativePosition(2)*combinedRadius;
                                -relativePosition(1)*combinedRadius + relativePosition(2)*leg]) / distSq;
                        end
                        
                        dotProduct2 = obj.dot(relativeVelocity,lineObj.direction);
                        u = dotProduct2*lineObj.direction - relativeVelocity;
                    end
                    % COLLISION
                else
                    % Collision. Project on cut-off circle of time timeStep.
                    invTimeStep = 1/dt;
                    % Vector from cutoff center to relative velocity. */
                    w = relativeVelocity - invTimeStep*relativePosition;
                    wLength = sqrt(obj.dot(w,w));
                    unitW = w/wLength;
                    
                    lineObj.direction = [unitW(2);-unitW(1)];
                    u = (combinedRadius*invTimeStep - wLength)*unitW;
                end
                % lineObj.point = va + 0.5*u;  % Equation 5 (definition of ORCA_AB)
                lineObj.point  = -va + 0.5*u;  % Equation 5 (definition of ORCA_AB)
                agnt_orcalines = [agnt_orcalines;lineObj];
            end
        end
        
        % LINEAR PROGRAM ONE (CORRECT)
        function [obj,flag,result_LP1] = linearProgram1(obj,lineVector,lineNo,constraintRadius,optVelocity,directionOpt,currentResult)
            % Solves a one-dimensional lineObjar program on a specified lineObj
            % subject to linear constraints defined by lineObjs and a circular
            % constraint.
            % INPUTS:
            % obj          - The agent object
            % lineObjs        - A vector of lineObjs defining the lineObjar constraints.
            % lineObjNo       - The specified lineObj constraint
            % radius       - The radius of the circulat constraint.
            % optVelocity  - The optimization velocity.
            % directionOpt - True if the direction should be optimized.
            % result       - A reference to the result of the lineObjar program.
            % OUTPUT:
            % flag - True if successful.
            
            
            % THE PROBLEM IS THAT THE INDEX NUMBER 'lineNo' is different
            % between linearProgram1 and linearProgram2
            
            
            dotProduct   = dot(lineVector(lineNo).point,lineVector(lineNo).direction);
            %            discriminant = dotProduct^2 + constraintRadius^2 - obj.absSq(lineVector(lineNo).point);
            discriminant = dotProduct^2 + constraintRadius^2 - dot(lineVector(lineNo).point,lineVector(lineNo).point);
            
            % Check if max speed circle fully invalidates lineObj constraint "lineObjNo"
            if discriminant < 0
                flag = logical(false);
                result_LP1 = currentResult;
                return
            end
            
            sqrtDiscriminant = sqrt(discriminant);
            tLeft  = -dotProduct - sqrtDiscriminant;
            tRight = -dotProduct + sqrtDiscriminant;
            
            % CHANGED TO PREVENT LOOP FROM ENTERING WHEN THE FAILURE OCCURS
            % ON THE FIRST ORCA-LINE
            i = 1;
            while (i < lineNo)
                denominator = obj.det(lineVector(lineNo).direction,lineVector(i).direction);
                numerator   = obj.det(lineVector(i).direction,(lineVector(lineNo).point - lineVector(i).point));
                
                if sqrt(dot(denominator,denominator)) <= obj.RVO_EPSILON    % The norm is less than the RVO constant
                    % lineObjs "lineObjNo" and "i" are (almost) parallel.
                    if numerator < 0
                        flag = logical(false);
                        result_LP1 = currentResult;
                        return
                    else
                        i = i + 1;                                          % Increment the counter before moving to next loop
                        continue
                    end
                end
                
                t = numerator/denominator;
                if denominator >= 0
                    tRight = min([tRight;t]);                               % lineObj "i" bounds lineObj "lineObjNo" on the right
                else
                    tLeft = max([tLeft,t]);                                 % lineObj "i" bounds lineObj "lineObjNo" on the left
                end
                
                if tLeft > tRight
                    flag = logical(false);
                    result_LP1 = currentResult;
                    return
                end
                i = i + 1;                                                  % Increment the counter before moving to next loop
            end
            
            % OPTIMISE DIRECTION & VELOCITY
            if directionOpt
                % Optimize direction
                if dot(optVelocity,lineVector(lineNo).direction) > 0
                    % Take right extreme
                    currentResult = lineVector(lineNo).point + tRight*lineVector(lineNo).direction;
                else
                    % Take left extreme
                    currentResult = lineVector(lineNo).point + tLeft*lineVector(lineNo).direction;
                end
            else
                % Optimise closest point
                t = dot(lineVector(lineNo).direction,(optVelocity - lineVector(lineNo).point));
                if t < tLeft
                    currentResult = lineVector(lineNo).point + tLeft*lineVector(lineNo).direction;
                elseif t > tRight
                    currentResult = lineVector(lineNo).point + tRight*lineVector(lineNo).direction;
                else
                    currentResult = lineVector(lineNo).point + t*lineVector(lineNo).direction;
                end
            end
            % DEFINE OUTPUT PARAMETERS
            result_LP1 = currentResult;
            flag = logical(true);
            
%             obj.LP1 = currentResult; % <<< REMOVE ME LATER
        end
        % LINEAR PROGRAM TWO
        function [obj,lineObjFlag,result_LP2]  = linearProgram2(obj,lineVector,constraintRadius,optVelocity,directionOpt)
            % Solves a two-dimensional lineObjar program subject to lineObjar
            % constraints defined by lineObjs and a circular constraint.
            % INPUTS:
            % obj          - The agent object
            % lineObjs        - lineObjs defining the lineObjar constraints.
            % radius       - The radius of the circular constraint.
            % optVelocity  - The optimization velocity.
            % directionOpt - True if the direction should be optimized.
            % result       - A reference to the result of the lineObjar program.
            % OUTPUT:
            % size_t       - The number of the lineObj it fails on, and the number of lineObjs if successful.
            
            if directionOpt
                % Optimise direction. Noe that the optimzation velocity is
                % of unit length in this case.
                currentResult = optVelocity*constraintRadius;
            elseif dot(optVelocity,optVelocity) > constraintRadius^2
                % Optimize closest point and ouside circle.
                currentResult = obj.normalise(optVelocity)*constraintRadius;
            else
                % Optimise closest point and inside circle
                currentResult = optVelocity;
            end

            % MOVE THROUGH THE lineObjS AND COMPARE THE lineObjAR CONSTRAINTS
            for i = 1:numel(lineVector)
                
                if obj.det(lineVector(i).direction,(lineVector(i).point - currentResult)) > 0
                    % "result" does not satisfy constraint "i". Compute new optimal "result".
                    tempResult = currentResult;
                    
                    % EVALUATE lineObj USING
                    [obj,flag,currentResult] = obj.linearProgram1(lineVector,i,constraintRadius,optVelocity,directionOpt,currentResult);
                    
                    if ~flag
                        % RETURN THE TEMP RESULT
                        result_LP2 = tempResult;
                        lineObjFlag = i;
                        return
                    end
                end
                
            end
            
            % RETURN THE FINAL RESULT OF THE LINEAR PROGRAM
            result_LP2 = currentResult;
            % RETURN THE TOTAL NUMBER OF LINES
            lineObjFlag = numel(lineVector); % The number of the line it fails on, and the number of lines if successful.
        end
        % LINEAR PROGRAM THREE
        function [obj,result_LP3] = linearProgram3(obj,lineVector,numObstLines,beginLine,constraintRadius,currentResult)
           % Solves a two-dimensional lineObjar program subject to lineObjar
           % constraints defined by lineObjs and a circular constraint.
           % INPUTS:
           % obj             - The agent object
           % lineObjs        - Linears defining the lineObjar constraints.
           % numObstlineObjs - Count of obstacle lineObjs.
           % beginlineObj    - The lineObj on which the 2-d linear program failed.
           % radius          - The radius of the circular constraint.
           % result          - A reference to the result of the linear program.
           
           % INITIAL DISTANCE TO "BEST CASE" POINT
           distance = 0;
           
           % IF THERE ARE NO OBSTACLES, ONLY AGENTS
           if numObstLines == 0
               numObstLines = 1; % Correct for c++ indexing
           end
           
           for i = beginLine:numel(lineVector)
               
               if obj.det(lineVector(i).direction,(lineVector(i).point - currentResult)) > distance
                   % Results does not satisfy constraint of lineObj 1
                   projlineObjs = [];
                   for j = numObstLines:i
                       % CREATE LINE OBJECT
                       newLine = struct('point',[],'direction',[]);
                       
                       determinant = obj.det(lineVector(i).direction,lineVector(j).direction);
                       if sqrt(dot(determinant,determinant)) <= obj.RVO_EPSILON
                           % lineObj "i" and lineObj "j" are parallel.
                           if dot(lineVector(i).direction,lineVector(j).direction) > 0
                               % lineObj "i" and lineObj "j" point in the same
                               % direction
                               continue;
                           else
                               % lineObj "i" and lineObj "j" point in opposite
                               % directions
                               newLine.point = 0.5*(lineVector(i).point + lineVector(j).point);
                           end
                       else
                           newLine.point = lineVector(i).point + (obj.det(lineVector(j).direction,...
                                          (lineVector(i).point - lineVector(j).point)) / determinant) * lineVector(i).direction;
                       end
                       
%                        newLine.direction = obj.normalise(lineVector(j).direction - lineVector(i).direction);
                       direction = lineVector(j).direction - lineVector(i).direction;
                       newLine.direction = direction/norm(direction);
                       projlineObjs = [projlineObjs;newLine];              % Append the new line to the projection lines
                   end
                   % HOLD ONTO THE CURRENT RESULT AS THE TEMP RESULT
                   
                   tempResult = currentResult;
                   
                   flowVector = [-lineVector(i).direction(2); lineVector(i).direction(1)]; % Vector in the best direction "going with the flow"
                   
                   % RE-EVALUATE THE NEW LINE OBJECTS
                   [obj,lineObjFlag,currentResult] = obj.linearProgram2(projlineObjs,constraintRadius,flowVector,logical(true));
%                    obj.LP2 = currentResult;
                   if lineObjFlag < numel(projlineObjs)
                       % This should in principle not happen.  The result is by definition
                       % already in the feasible region of this linear program. If it fails,
                       % it is due to small floating point error, and the current result is
                       % kept.
                       currentResult = tempResult;
                   end
                   % The distance to the new optimal flow line
                   distance = obj.det(lineVector(i).direction,(lineVector(i).point - currentResult));
               end
           end
           % REDEFINE THE RESULT
           result_LP3 = currentResult;
%            obj.LP3 = result_LP3;
       end
    end
    %/////////////////////// "agent.cpp" FUNCTIONS ////////////////////////
    methods (Static)
        % AGENT UPDATE FUNCTION (CORRECT)
        function [position,velocity] = update(position,newVelocity)
            % INPUT HANDLING
            assert(numel(position) == numel(newVelocity),'Position and velocity updates must be of the same dimensions.');
            % UPDATE THE AGENT'S VELOCITY
            velocity = newVelocity;
            position = position + dt*velocity; % Update in accordance to the RVO2 library
        end
    end
    
    % ////////////////////// "vector_2.h" FUNCTIONS ///////////////////////
    methods (Static)
        % COMPUTES THE LENGTH OF A 2D VECTOR
        function [absLength] = abs(v)
           absLength = sqrt(dot(v,v));   
        end
        % RETURN THE SQUARE LENGTH OF A 2D VECTOR
        function [absSq] = absSq(v)
           absSq = dot(v,v);
        end
        % COMPUTE THE DETERMININANT OF A 2D SQUARE MATRIX
        function [det_uv] = det(u,v)
            det_uv = u(1)*v(2) - u(2)*v(1);
        end
        % THE DOT PRODUCT
        function [dot_uv] = dot(u,v)
           dot_uv = u(1)*v(1) + u(2)*v(2); 
        end
        % NORMALIZE A VECTOR
        function [unitVector] = normalise(v)
           unitVector = v/sqrt(dot(v,v));
        end
    end
    
    %/////////////////////// HELPER FUNCTIONS [JD] ////////////////////////
    methods (Static)
        % BUILD THE VERTEX-OBSTACLE LIST
        function [vertexSet] = defineRVOVertexObstacleSet(vertexData)
            % This function takes a list of points and contstucts a linked
            % obstacle list description. The vertices are assumed ordered
            % in terms of their connectivity (i.e. clockwise).
            
            
            % VERTEX STRUCTURE
            obstacleTemp = struct(...
                'point',zeros(2,1),...
                'nextObstacle',0,...
                'prevObstacle',0,...
                'unitDir',zeros(2,1),...
                'isConvex',boolean(true),...
                'id',0);
            
            numberOfVertices = size(vertexData,1);
            
            % LIST OF OBSTACLES
            vertexSet = repmat(obstacleTemp,numberOfVertices, 1 );
            % MOVE THROUGH THE VERTEX SET
            for i = 1:numberOfVertices
                % INDEX OF i+1
                ind_i_plus = mod(i,numberOfVertices) + 1;
                ind_i_plus2 = mod(i+1,numberOfVertices) + 1;
                % ALLOCATE POINT
                vertexSet(i).point = vertexData(i,:)';                     % Assign obstacle
                % SET NEXT AND PREVIOUS OBSTACLES/VERTEX POINTS
                vertexSet(i).nextObstacle = ind_i_plus;
                vertexSet(ind_i_plus).prevObstacle = i;
                % UNIT DIRECTION
                lineDirection = vertexData(ind_i_plus,:)' - vertexData(i,:)';
                vertexSet(i).unitDir = OMAS_geometry.unit(lineDirection);
                % IS CONVEX
                if numberOfVertices == 2
                    vertexSet(ind_i_plus).isConvex = 1;
                else
                    % CHECK IF NEXT POINT IS ON LEFT (IS CONVEX)
                    isLeft = OMAS_geometry.leftOf(vertexData(i,:),vertexData(ind_i_plus,:),vertexData(ind_i_plus2,:));
                    % IS POINT LEFT OF
                    vertexSet(ind_i_plus).isConvex = isLeft <= 0;          % If list is clockwise "<=0" , else list is not ">=0"
                end
                % SET THE INDEX OF POINT
                vertexSet(i).id = i;
            end      
        end
    end
end