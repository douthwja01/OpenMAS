%% THE 2D AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the 'agent' class, this object class is intended provide
% support for 2D simulation. The functions provided here 

% Author: James A. Douthwaite

classdef agent_2D < agent
%%% AGENT(2D) BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties

    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_2D(varargin)
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
            obj@agent(varargin); % Call the super class
            % VIRTUAL DEFINITION
            obj.VIRTUAL.type = OMAS_objectType.agent;
            obj.VIRTUAL.symbol = 'diamond';
            obj.VIRTUAL.detectionRange = obj.sensorRange;                  % Assume the agent has perfect environmental knowledge (m)
            % AGENT DEFAULT KINEMATIC LIMITATIONS
            obj.linearVelocityLimits = [inf;inf;inf];
            obj.angularVelocityLimits = [inf;inf;inf];
            obj.linearAccelerationLimits = [inf;inf;inf];
            obj.angularAccelerationLimits = [inf;inf;inf];                 % Emulating no limitations
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
            else
                error('Object TIME packet is invalid.');
            end
            
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            observationSet = varargin{1}; % The detected objects
            [obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(observationSet);

            % GET WAYPOINT (TARGET) HEADING VECTOR \\\\\\\\\\\\\\\\\\\\\\\\
            if ~isempty(obj.targetWaypoint)  
               % GET THE WAYPOINT HEADING VECTOR
               targetHeading = obj.targetWaypoint.state(1:2)/norm(obj.targetWaypoint.state(1:2));
            else
               targetHeading = [1;0];
            end
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            linearAcceleration = [0.1;0];
            angularAcceleration = 0;
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_accelerations(dt,linearAcceleration,angularAcceleration);
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
        % /////////////////////////////////////////////////////////////////
        
        % //////////////// GENERAL INTERACTIONS ///////////////////////////
        % SIMPLE CONTROLLER
        function [d_heading,d_speed,obj] = simpleController(obj,targetHeading,targetSpeed)
            % This function assumes that there is no controller and simply
            % computes the change in the agents current heading and speed
            % to meet the target reference.
            % INPUTS:
            % targetHeading - Target heading vector
            % targetSpeed   - Target scalar speed
            
            currentVelocity = obj.localState(7:8);                         % Agent's current velocity
            targetHeading = targetHeading/norm(targetHeading);             % Renormalise the heading
            [~,e_lambda] = obj.getControlInputs([1;0],targetHeading);      % Get the equivalent angular heading and elevation
            d_heading = -e_lambda;                                         % Instantaneous change in heading
            d_speed = targetSpeed - norm(currentVelocity);                 % Speed change is the speed difference
        end
        % IMPOSE 2D KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = kinematicContraints(obj,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            % HEADING IS RELATIVE
            currentVelocity = obj.localState(7:8);
            currentHeading = zeros(size(desiredHeading));
            
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - currentHeading;
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = obj.boundValue(Hdot,-obj.angularVelocityLimits,obj.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = currentHeading + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed;                          % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = obj.boundValue(acceleration,-obj.linearAccelerationLimits(1),obj.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = obj.boundValue(achievedSpeed,-obj.linearVelocityLimits(1),obj.linearVelocityLimits(1)); 
        end
        
        % ///////////// LOCAL(AGENT)(2D) CORE OPERATIONS //////////////////
        % ALL IN ONE AGENT-INFORMATION UPDATE (FOR 2D STATES)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment in 2D.
            % INPUTS:
            % dt              - The unit timstep
            % observedObjects - The full observed object structure
            % OUTPUTS:
            % obj             - The updated object
            % obstacleSet     - The list of obstacle structures
            % waypointSet     - The list of waypoint structures
            
            agentSet = [];
            obstacleSet = [];
            waypointSet = [];   % No observed objects 
            
            % INPUT HANDLING
            if isempty(observedObjects)
                return
            end
            positionIndices = 1:2;
            
            for item = 1:length(observedObjects)
                % DEFINE THE SENSED VARIABLES
                objectName = observedObjects(item).name;
                objectID   = observedObjects(item).objectID;               % The objects ID
                objectType = observedObjects(item).type;                   % The objects type (obstacle or waypoint)
                objectSize = observedObjects(item).size;
                objectState = [observedObjects(item).position(positionIndices,1);...
                               observedObjects(item).velocity(positionIndices,1)];     % Obstacle state [x y u v]
                objectTTC = observedObjects(item).timeToCollision;
                objectPriority = observedObjects(item).priority;
                
                % DEFINE ITS APPARENT PRIORITY
                if ~isnan(objectPriority)
                    objectPriority = observedObjects(item).priority;           % Assign priority to memory
                else
                    objectPriority = 1/objectTTC; 
                end
                
                % MEMORY ITEM CONSTRUCTIONS
                memoryItem = struct('name',objectName,'objectID',objectID,...
                                    'type',objectType,'size',objectSize,...
                                    'state',objectState,'TTC',observedObjects(item).timeToCollision,...
                                    'priority',objectPriority); 
                
                % UPDATE THE AGENTS KNOWLEDGE
                obj = obj.updateAgentKnowledge(memoryItem);
            end
            
            % SORT OBJECT SET BY PRIORITY
            [~,ind] = sort([obj.memory.priority],2,'descend');             % Ordered indices of the object IDs
            obj.memory = obj.memory(ind);
            
            % DECERN OBJECT TYPES
            % obj.memory now contains an updated understanding of all the
            % objects in its visual horizon. From this, the waypoints and
            % obstacles must be parsed.
            agentIndices    = [obj.memory.type] == OMAS_objectType.agent;
            waypointIndices = [obj.memory.type] == OMAS_objectType.waypoint;  
            obstacleIndices = [obj.memory.type] == OMAS_objectType.obstacle; % Differenciate between the different object types
            
            % PULL MEMORY ITEMS FOR LATER MANIPULATION
            agentSet    = obj.memory(agentIndices);
            waypointSet = obj.memory(waypointIndices);
            obstacleSet = obj.memory(obstacleIndices); 
            
            % UPDATE THE TARGET WAYPOINT, PRIORTIES 
            [obj,~] = obj.waypointUpdater(waypointSet);                    % Update the target waypoint and heading vector
        end
        % SELECT WAYPOINT (FOR 2D STATES)
        function [obj,waypointVector] = waypointUpdater(obj,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if ~exist('waypointSet','var') || isempty(waypointSet)
                obj.targetWaypoint = [];       % No waypoints are available
                waypointVector = [];           % Default heading
                return
            end
            
            positionIndices = 1:2;
            
            % WAYPOINT SET IS POPULATED -> GET THE PRIORITY VECTOR
            priorityVector = [waypointSet.priority];                       % The vector of priorities
            waypointIndex  = linspace(1,length(priorityVector),length(priorityVector)); % Memory of their position in the waypoint set
            waypointIDset  = [waypointSet.objectID];                       % The vector of IDs 
            waypointMatrix = [waypointIndex;waypointIDset;priorityVector];
            
            % WAYPOINT MATRIX IS OF THE FORM
            % [ 1 2 3 ] - Positions in the waypoint set.
            % [ 2 7 4 ] - Waypoint IDs.
            % [ 3 5 9 ] - Waypoint priorities.
                        
            % TARGET WAYPOINT IS EMPTY -> POPULATE WITH NEW WAYPOINTS
            if isempty(obj.targetWaypoint)
                % NO WAYPOINTS HAVE BEEN SELECTED
                if isempty(obj.achievedWaypoints)
                    [~,maxPriorityIndex] = max(waypointMatrix(3,:));       % Index of maximum priority value in waypoint matrix
                    waypointSetIndex = waypointMatrix(1,maxPriorityIndex); % Gets the associated waypoint set index of the highest priority
                    obj.targetWaypoint = waypointSet(waypointSetIndex);    % Update with the highest priority waypoint
                else
                    % SELECT HIGHEST PRIORITY, BUT NOT ACHIEVED
                    invalidWaypointIDs = obj.achievedWaypoints;            % Get the achieved waypoint IDs
                    availableWaypointIDs = waypointMatrix(2,:);            % Get those visable IDs
                    validWaypoints = availableWaypointIDs ~= invalidWaypointIDs;   % Get those visible IDs that are valid
                    % IF NO FURTHER VALID IDs
                    if ~any(validWaypoints)
                        obj.targetWaypoint = [];
                        waypointVector = [];
                        return
                    else
                        % SELECT REMAINING WAYPOINTS
                        validWaypoints = waypointMatrix(:,validWaypoints); % Get waypoinSet location, IDs and priority of valid waypoints
                        % SELECT WAYPOINT OF NEXT HIGHEST PRIORITY                        
                        [~,maxValidIndex ] = max(validWaypoints(3,:));       % Get the index of max priority waypoint
                        waypointSetIndex = validWaypoints(1,maxValidIndex);  % Get the location of the target waypoint in the waypoint set
                        obj.targetWaypoint = waypointSet(waypointSetIndex);  % Select waypoint object
                    end
                end
            end
            
            % EVALUATE CURRENT TARGET WAYPOINT ////////////////////////////          
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            waypointTolerance = 0.1;
            waypointAchievedCondition = norm(obj.targetWaypoint.state(positionIndices)) - (obj.targetWaypoint.size + obj.VIRTUAL.size);
            % IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
            if waypointAchievedCondition <= waypointTolerance && ~any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID))
                obj.achievedWaypoints = horzcat(obj.achievedWaypoints,obj.targetWaypoint.objectID); % Add achieved objectID to set
            end
            
            % CHECK IF CURRENT TARGET IS VALID ////////////////////////////
            % ASSESS THE ACHIEVED WAYPOINT SET AGAINST THE CURRENT TARGET
            invalidTargetCondition = any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID)); % Is the current target ID in the achieved set
            if invalidTargetCondition
               % TARGET WAYPOINT HAS BEEN ACHIEVED, REALLOCATE ONE
               selectionVector = priorityVector < obj.targetWaypoint.priority;   % Only those with lower priority
               if ~any(selectionVector)
                   obj.targetWaypoint = [];
                   waypointVector = [];
                   return
               else
                   reducedMatrix = waypointMatrix(:,selectionVector);            % Build vector of priorities less than 
                   [~,reducedIndex] = max(reducedMatrix(3,:));                   % Finds the max priority index where less than current
                   waypointIndex = reducedMatrix(1,reducedIndex);                % Get the waypoint index from the reduced priority matrix               
                   obj.targetWaypoint = waypointSet(waypointIndex);              % Select the next highest priority waypoint
               end
            else          
                currentWaypointID = obj.targetWaypoint.objectID;                 % Get the ID of the current target waypoint
                selector = (waypointMatrix(2,:) == currentWaypointID);           % Find where the ID appears in the waypoint matrix
                waypointIndex = waypointMatrix(1,selector);                      % Get the waypoint Set index from the waypoint matrix
                obj.targetWaypoint = waypointSet(waypointIndex);                 % Update with the highest priority waypoint  
            end           
            % THE CORRECT TARGET WAYPOINT IS NOW DEFINED, GENERATE HEADING                      
            waypointVector = obj.targetWaypoint.state(positionIndices)/norm(obj.targetWaypoint.state(positionIndices));  
        end
        % STATE UPDATES & DYNAMICS
        % The basic utilities for updating the local state vector are
        % handled in the agent class and operate in both the 2D and 3D
        % case.
        
    end
    %% STATIC AGENT TOOLS
    methods (Static)
        % DEFINE THE 2D CONTROL INPUTS FROM DESIRED (2D) VELOCITY
        function [theta,lambda] = getControlInputs(V,U)
            % INPUTS
            % V - Current velocity vector [u v] (comparative vector)
            % U - The unit correction vector [u v]
                        
            % NORMALISE THE ELEMENTS
            V = V/norm(V); U = U/norm(U); % Normalise vector 
            % GET THE HORIZONTAL ELEMENTS
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];                                            % Reject the vertical elements
            % DEFINE LOS ANGLE            
            rotationAxis = cross(Vh,Uh);
            lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            % GET THE ELEVATION ANGLE 
            theta = atan2(Uh(3),norm(Uh));
        end
        % CALCULATE THE NEW (2D) STATE ESTIMATE
        function [newState] = linearStateEstimation(priorState,newPosition,dt)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
           % INPUTS:
           % priorState  - [Sx;Sy;Vx;Vy];
           % newPosition - The new relative position
           % OUTPUT:
           % newState - The new obstacle state
 
           % GET THE PRIOR VALUES
           position0 = priorState(1:2);
           velocity0 = priorState(3:4);
           % GET THE POSITION
           position1 = newPosition;
           dX = (position1 - position0);
           velocity1 = dX/dt;  % Defines the average velocity
           % NO PREVIOUS VELOCITY RECORDED
           if any(isnan(velocity0))
              newState = [position1;velocity1];
              return
           end
           
           % ESTIMATE THE CHANGE IN VELOCITY
           newState = [newPosition;velocity1];  % Append to the obstacles state
        end
    end  
end

