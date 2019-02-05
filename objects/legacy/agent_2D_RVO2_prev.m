%% THE RECIPROCAL VELOCITY OBSTACLE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This agent class provides support for the RVO2 library for ORCA-based
% Reciprocal Velocity Obstacle collision avoidance.

% Author: James A. Douthwaite

classdef agent_2D_RVO2 < agent_2D
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % AGENT PARAMETERS
        maxSpeed  = 2;       % Defined in the example defaults
        nominalSpeed = 1;
        radius = 0.5;        % Declare (independant of the VO varients)
        % AVOIDANCE PARMETERS
        RVO_EPSILON   = 1E-5;
        neighbourDist = 15;  %(s)
        maxNeighbours = 10;
        timeHorizon   = 10;
        timeHorizonObst = 10;
        orcaLines = struct('point',[],'direction',[]);
        
        % DIAGNOSITICS
        TIME = [];
        testVariable = NaN(2,1);
        LP1 = NaN(2,1);
        LP2 = NaN(2,1);
        LP3 = NaN(2,1);
        
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
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D(varargin); 
            
            % GET THE SENSOR
            obj = obj.getSensorParameters();
            
            % VIRTUAL DEFINITION (SIMULATOR PARAMETERS)
            obj.VIRTUAL.radius = obj.radius; 
            obj.VIRTUAL.detectionRange = inf;
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        % ///////////////// AGENT MAIN CYCLE //////////////////////////////
        function obj = processTimeCycle(obj,TIME,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
%                 obj.TIME = TIME; %    <<< REMOVE ME
            else
                error('Object TIME packet is invalid.');
            end
            
            % DEFAULT BEHAVIOUR 
            desiredSpeed = obj.nominalSpeed;
            desiredHeadingVector = [1;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
                        
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});       % IDEAL INFORMATION UPDATE
%             [obj,obstacleSet,agentSet] = obj.getSensorUpdate(dt,varargin{1}); % REALISTIC INFORMATION UPDATE

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
                desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
            end
                       
            % ////////// COLLISION AVOIDANCE VECTOR GENERATION ////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(obstacleSet) || ~isempty(agentSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE AGENTS UPDATED PREFERRED VELOCITY
                [obj,desiredVelocity] = obj.getAvoidanceCorrection(dt,desiredVelocity,agentSet,obstacleSet); 
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
%             % CONVERT THE DESIRED VELOCITY TO HEADING
%             [lambda] = obj.getVectorHeadingAngles([1;0;0],[desiredVelocity;0]);
%             newHeadingRate = -lambda/dt; % -ve feedback
%             
%             obj.testVariable = newHeadingRate;
            
            % APPLY SPEED CONSTRAINT
            if norm(desiredVelocity) > obj.maxSpeed
                desiredVelocity = (desiredVelocity/norm(desiredVelocity))*obj.maxSpeed;
            end    
            
            % /////////////////// STATE DYNAMICS //////////////////////////
%             [newState] = obj.stateDynamics_simple(dt,[norm(desiredVelocity);0],newHeadingRate);
%             obj = obj.updateGlobalProperties_fixedFrame(dt,newState); 
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////   
%             obj = obj.updateGlobalProperties(dt,newState); 
             
            % FOR EVALULATION:
            [newState] = obj.update(dt,desiredVelocity);    
            newState = [newState(1);newState(2);0;newState(3);newState(4);0]; % Map to 2D state vector
            
            obj = obj.updateGlobalProperties_fixedFrame(dt,newState);  

            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(TIME,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'Vx (m/s)','Vy (rad)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [newState(4);newState(5);newState(6)];         % Record the control inputs    
        end
    end
    % //////////////////////// OPENMAS FUNCTIONS ///////////////////////////
    methods 
        % SPECIFIC RVO STATE INITIALISER
        function obj = initialise_localState(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y phi dx dy dphi]
            
            % Build the initial state vec
            obj.localState = zeros(6,1);
            obj.localState(3) = localXYZrotations(3);     % The initial heading angle
            obj.localState(4:5) = localXYZVelocity(1:2);  % The initial 2D velocities  
            % ADDITIONAL PARAMETERS
%             obj.position = obj.localState(1:2,1);
%             obj.velocity = obj.localState(4:5,1);
        end
        % MAPPING OF OPENMAS TO THE RVO LIBRARY
        function [obj,avoidanceVelocity] = getAvoidanceCorrection(obj,dt,prefVelocity,agentSet,obstacleSet)
            % This function computes the mapping of the simulation inputs
            % for the RVO library functions and returns the calculated
            % optimal velocity.
            
            % ADDITIONAL PARAMETERS
            agentPosition = obj.localState(1:2,1);
            agentVelocity = obj.localState(4:5,1);
            
            % BUILD THE REQUIRED AGENT AND OBSTACLE STRUCTURES
            agentNeighbours = [];
            agentNeighbourLength = 0;
            for item = 1:numel(agentSet)
                % NEIGHBOUR CONDITIONS
%                 neighbourConditionA = agentNeighbourLength < obj.maxNeighbours;            % Maximum number of neighbours
%                 neighbourConditionB = norm(agentSet(item).position) < obj.neighbourDist;      % [CONFIRMED] 
%                 if ~neighbourConditionA || ~neighbourConditionB
%                     continue
%                 end

                % NEIGHBOUR CONDITIONS
                neighbourConditionB = norm(agentSet(item).position) < obj.neighbourDist;  % [CONFIRMED] 
                neighbourConditionC = ~any(isnan(agentSet(item).velocity));               % Wait for a valid velocity reading
                if ~neighbourConditionB || ~neighbourConditionC
                    continue
                end
                
                % BUILD THE RVO AGENT DESCRIPTION
                neighbourPosition = agentSet(item).position + agentPosition;
                neighbourVelocity = agentSet(item).velocity + agentVelocity;
                agentRef = struct('position',neighbourPosition,...
                                  'velocity',neighbourVelocity,...
                                  'unitDir',neighbourVelocity/norm(neighbourVelocity),...
                                  'radius',agentSet(item).radius,...
                                  'id',agentSet(item).objectID);
                agentNeighbours = [agentNeighbours;agentRef];
                agentNeighbourLength = numel(agentNeighbours);
            end
            
            obstacleNeighbours = [];
            obstacleNum = numel(obstacleSet);
            for item = 1:obstacleNum
                k = mod(obstacleNum,item);
                if k == 0 && item == 1 % BEGINNING OF LIST
                    nextObstacle = obstacleSet(item+1);
                    prevObstacle = obstacleSet(obstacleNum);
                elseif k == 0 && item ~= 1   % END OF LIST
                    nextObstacle = obstacleSet(1);
                    prevObstacle = obstacleSet(item-1);
                else
                    nextObstacle = obstacleSet(item+1);
                    prevObstacle = obstacleSet(item-1);
                end
                % BUILD NEIGHBOUR OBJECT
                nextNeighbour = struct('isConvex',logical(false),...
                                  'nextObstacle',[],...
                                  'prevObstacle',[],...
                                  'point',nextObstacle.state(1:2,1),...
                                  'unitDir',nextObstacle.state(3:4,1)/norm(nextObstacle.state(3:4,1)),...
                                  'id',nextObstacle.objectID);
                prevNeighbour = struct('isConvex',logical(false),...
                                  'nextObstacle',[],...
                                  'prevObstacle',[],...
                                  'point',prevObstacle.state(1:2,1),...
                                  'unitDir',prevObstacle.state(3:4,1)/norm(prevObstacle.state(3:4,1)),...
                                  'id',prevObstacle.objectID);
                
                obstacle1 = obstacleSet(item);
                obstacleRef = struct('isConvex',logical(false),...
                                  'nextObstacle',nextNeighbour,...
                                  'prevObstacle',prevNeighbour,...
                                  'point',obstacleSet(item).state(1:2,1),...
                                  'unitDir',obstacle1.state(3:4,1)/norm(obstacle1.state(3:4,1)),...
                                  'id',obstacle1.objectID);                % As defined in obstacle.h
                obstacleNeighbours = [obstacleNeighbours;obstacleRef];
            end
                        
            % COMPUTE THE AGENT'S NEW VELOCITY
            % This calls the libraries 'computeNewVelocity' function
            [obj,avoidanceVelocity] = obj.computeNewVelocity(dt,prefVelocity,agentNeighbours,obstacleNeighbours);
        end
        % INITIALISE SENSORS (NOISY SENSOR PARAMETERS)
        function [obj] = getSensorParameters(obj)
            % This function is designed to populate the SENSOR field with
            % representative sensor uncertainty.
            % BUILD THE SENSOR DEFINITION
            obj.SENSORS = struct('sensorHorizon',inf,...    % Assume the agent has perfect environmental knowledge (m)
                                 'positionSigma',0.5,...    % Accurate to within 0.5m
                                 'velocitySigma',0.1,...    % Accurate to within 0.1m/s
                              'rangeFinderSigma',0.1,...    % Accurate to within 0.1m
                                   'cameraSigma',5.21E-5,...% One pixel in a 1080p image
                               'sampleFrequency',inf);      % Sensing has perfect precision
        end
    end
    %% FUNCTIONS FROM THE "agent.cpp" FILE
    methods
       % COMPUTE THE AGENT VELOCITY
       function [obj,avoidanceVelocity] = computeNewVelocity(obj,dt,prefVelocity,agentNeighbours,obstacleNeighbours)
            % This function computes the new velocity of this agent from the
            % RVO2 library.
            
            % DEFINE THE AGENT PARAMETERS
            agentPosition = obj.localState(1:2,1);
            agentVelocity = obj.localState(4:5,1);
            agentRadius   = obj.VIRTUAL.radius;
            
            % MAP THE PREFERRED VELOCITY TO NEW FRAME
            prefVelocity = -prefVelocity;
            
            % RESET THE OCRA-LINES VECTOR FROM THE LAST TIMESTEP
            obj.orcaLines = [];
            
            % //////////// GET THE OBSTACLE ORCA LINES ////////////////////
            invTimeHorizonObst = 1/obj.timeHorizonObst;
            if ~isempty(obstacleNeighbours)
            end
            numObstLines = numel(obj.orcaLines);
            
            % ////////////// GET THE AGENT ORCA LINES /////////////////////
            % BEGIN CREATING ORCA LINES FOR THE AGENTS
            invTimeHorizon  = 1/obj.timeHorizon;

            for i = 1:numel(agentNeighbours)
                % CALCULATE THE OTHER AGENT PROPERTIES
                other = agentNeighbours(i);
%                     relativePosition = other.position - agentPosition;
%                     relativeVelocity = agentVelocity - other.velocity;
                
                relativePosition = -(other.position - agentPosition);
                relativeVelocity = -(agentVelocity - other.velocity);
                
                distSq = dot(relativePosition,relativePosition);
                combinedRadius = agentRadius + other.radius;
                combinedRadiusSq = combinedRadius^2;
                
                % CREATE ORCA lineObj
                lineObj = struct('point',[],'direction',[]);
                
                % NO COLLISION 
                if distSq > combinedRadiusSq
                    w = relativeVelocity - invTimeHorizon*relativePosition; % Vector from cutoff center to relative velocity
                    wLengthSq   = dot(w,w); % Square length
                    dotProduct1 = dot(w,relativePosition);
                    
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
                        
                        dotProduct2 = dot(relativeVelocity,lineObj.direction);
                        u = dotProduct2*lineObj.direction - relativeVelocity;
                    end
                % COLLISION     
                else
                    % Collision. Project on cut-off circle of time timeStep.
                    invTimeStep = 1/dt;
                    % Vector from cutoff center to relative velocity. */
                    w = relativeVelocity - invTimeStep*relativePosition;
                    wLength = sqrt(dot(w,w));
                    unitW = w/wLength;
                    
                    lineObj.direction = [unitW(2);-unitW(1)];
                    u = (combinedRadius*invTimeStep - wLength)*unitW;
                end
%                 lineObj.point = agentVelocity + 0.5*u;  % Equation 5 (definition of ORCA_AB)
                lineObj.point = -agentVelocity + 0.5*u;  % Equation 5 (definition of ORCA_AB)                
                obj.orcaLines = [obj.orcaLines;lineObj];
            end
            
            % RESET THE DIAGNOSTIC VARIABLES
%             if obj.objectID == 1
%                 obj.LP1 = NaN(2,1);
%                 obj.LP2 = NaN(2,1);
%                 obj.LP3 = NaN(2,1);
%                 obj.testVariable = NaN(2,1);
%             end
            
%             if obj.TIME.currentTime == 9.5
%                display('problem time'); 
%                obj.orcaLines
%             end
            
            % CALL THE SECOND LINEAR PROGRAM
            [obj,lineObjFail,optimalVelocity] = obj.linearProgram2(obj.orcaLines,obj.maxSpeed,prefVelocity,logical(false));
            
            % IF SOME ORCA-LINES FAIL
            % This function will generate direction based ORCA lines where
            % the previous lines failed.
            if lineObjFail < numel(obj.orcaLines) 
                [obj,optimalVelocity] = obj.linearProgram3(obj.orcaLines,numObstLines,lineObjFail,obj.maxSpeed,optimalVelocity);
            end 
            % REASSIGN OUTPUT VELOCITY                               
%             avoidanceVelocity = optimalVelocity;             
            avoidanceVelocity = -optimalVelocity;      
            
            % ////////// WRITE DATA TO CSV //////////////////////////
%             if obj.objectID == 1 
%                 params = [obj.TIME.currentTime,obj.objectID,...
%                           obj.VIRTUAL.globalPosition(1),...
%                           obj.VIRTUAL.globalPosition(2),...
%                           obj.VIRTUAL.globalVelocity(1),...
%                           obj.VIRTUAL.globalVelocity(2),...
%                           8888,...
%                           1000,obj.LP1(1),obj.LP1(2),...
%                           2000,obj.LP2(1),obj.LP2(2),...
%                           3000,obj.LP3(1),obj.LP3(2)];
% %                           8888,desiredVelocity(1),desiredVelocity(2),...
% %                           8888,obj.testVariable(1)];
%                 % WRITE TO CSV FILE
%                 obj.sendVectorToCsv(obj.TIME,params)
%             end         
       end
       % AGENT UPDATE FUNCTION (CORRECT)
       function [state,position,velocity] = update(obj,dt,newVelocity)
           % INPUT HANDLING
           assert(size(newVelocity,1) == 2*size(newVelocity,2),'New velocity has wrong dimensions');
           % EVALUATE STATE UPDATES
           velocity(1:2,1) = newVelocity;
           position(1:2,1) = obj.localState(1:2,1) + dt*velocity;
           state = [position;velocity];
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
           
           obj.LP1 = currentResult; % <<< REMOVE ME LATER
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
                      
           obj.LP2 = currentResult; % <<< REMOVE ME LATER      
           
%            if obj.TIME.currentTime == 9.5
%                display('problem time');
%                obj.orcaLines
%            end
           
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
                   obj.LP2 = currentResult;
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
           obj.LP3 = result_LP3;
       end
    end
    
    methods (Static)
        %% FUNCTIONS IN THE VECTOR_2.h FILE  
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
end
