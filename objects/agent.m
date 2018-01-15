%% THE AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic agent and import this variables 
% into the simulation space for the purpose of multi-vehicle control simulation.
% The agent object is a child of the objectDefintion; the prinicple
% distinctions being:
% sensorRange      - The agent is assumed capable of observing its
%                    surroundings.
% controlFrequency - The frequency at which the control cycle is computed.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef agent < objectDefinition
%%% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % PLACE-HOLDER KINEMATIC LIMITATIONS
        linearVelocityLimits;           % Limits on the agents linear velocity 
        linearAccelerationLimits;       % Limits on the agents linear acceleration
        angularVelocityLimits;          % Limits on the agents angular velocity
        angularAccelerationLimits;      % Limits on the agents angular acceleration
        % CONTROLLER PARAMETERS
        priorError = [];            % Container for prior control error
        % OBSTACLE OBJECT MATRIX (record of sightings)
        memory;                     % Record of last known [objectID;positions;velocities]
        % WAYPOINTS
        targetWaypoint;             % The current waypoint target
        achievedWaypoints;          % The agents list of locally achieved waypoints
        % SIMULATION IDENTIFIERS
        sensorRange = inf;          % Assume the agent has perfect environmental knowledge (m)
        % SAMPLE FREQUENCY
        sampleFrequency = inf;      % Object has perfect precision
        % AGENT-SIDE OUTPUT DATA
        DATA;                       % The output container for agent-side data.
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent(varargin)
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
            obj@objectDefinition(varargin); % Call the super class
            
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
        % AGENT MAIN CYCLE 
        function obj = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
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
            
            % CHECK FOR NEW INFORMATION UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\\\\
            observationSet = varargin{1}; % The detected objects
            
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(observationSet);

            % GET WAYPOINT (TARGET) HEADING VECTOR \\\\\\\\\\\\\\\\\\\\\\\\
            if ~isempty(obj.targetWaypoint)  
               % GET THE WAYPOINT HEADING VECTOR
               targetHeading = obj.targetWaypoint.state(1:3)/norm(obj.targetWaypoint.state(1:3));
            else
               targetHeading = [1;0;0];
            end
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            linearAcceleration = [0;0;0];
            angularAcceleration = [0;0;0];
            
            % GET THE NEW STATE VECTOR
            newState = obj.stateDynamics_accelerations(dt,linearAcceleration,angularAcceleration);
            % UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,newState);
        end
        
        % //////////////// GENERAL INTERACTIONS ///////////////////////////
        % FREEZE AGENT
        function [newState] = freezeAgent(obj)
            % This function zeros the agents linear and angular velocities
            newState = obj.localState;
%             newState(7:12) = zeros(6,1);
            newState(7:12) = NaN([6,1]);
        end  
        
        % /////////////// BASIC CONTROL MECHANISM /////////////////////////
        % SPEED AND HEADING (4D) PID CONTROLLER
        function [d_heading,d_speed,obj] = trajectoryController(obj,targetHeading,targetSpeed)
            % This function is desiged to generate feedack from a desired
            % local velocity vector. 
            % INPUTS:
            % targetHeading - The unit heading vector
            % targetSpeed   - The target speed in that direction
            % OUTPUTS:
            % heading_fb    - The feedback on the agent heading
            % speed_fb      - The feedback on the agent speed
            % obj           - The error-updated agent object
            
            % RENORMALISE THE HEADING VECTOR
            targetHeading = targetHeading/norm(targetHeading);
            
            % DEFINE THE STATE ERROR TERMS
            % RELATIVE HEADING
            [e_theta,e_lambda] = obj.getControlInputs([1;0;0],targetHeading);  % Get the equivalent angular heading and elevation
            % RELATIVE SPEED
            e_speed = targetSpeed - norm(obj.localState(7:9));                 % The speed error 
            controlError = [e_speed;0;e_theta;e_lambda];                       % -ve in the NED control frame of reference
            
            % TUNING PARAMETERS 
            Kp_linear = 0.8;
            Kd_linear = 0;
            Kp_angular = 1;
            Kd_angular = 0;
            
            % CALCULATE THE CONTROL FEEDBACK
            control_fb = diag([Kp_linear Kp_angular Kp_angular Kp_angular])*controlError + ...
                         diag([Kd_linear Kd_angular Kd_angular Kd_angular])*(controlError - obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = controlError;
            
            % PARSE THE CONTROL FEEDBACK FOR EXTERNAL USE
            d_speed    = control_fb(1);             % Absolute speed input
            d_heading  = control_fb(2:4);
            % CONVERT FLU (AVOIDANCE TO NED (CONTROL)
            d_heading = diag([1 1 -1])*d_heading;
        end
        % SIMPLE CONTROLLER
        function [d_heading,d_speed,obj] = simpleController(obj,targetHeading,targetSpeed)
            % This function assumes that there is no controller and simply
            % computes the change in the agents current heading and speed
            % to meet the target reference.
            % INPUTS:
            % 
            
            % RENORMALISE THE HEADING VECTOR
            targetHeading = targetHeading/norm(targetHeading);             % Re-normalise the target heading
            % RELATIVE HEADING
            [e_theta,e_lambda] = obj.getControlInputs([1;0;0],targetHeading);  % Get the equivalent angular heading and elevation
            d_heading = [0;e_theta;-e_lambda];                                 % Instantaneous change in heading
            % GET THE CHANGE TO BRING ABOUT DESIRED SPEED
            d_speed = targetSpeed - norm(obj.localState(7:9));                 % Speed change is the speed difference
        end
        % IMPOSE KINEMATIC LIMITS ON CONTROL INPUTS
        function [achievedHeading,achievedSpeed] = kinematicContraints(obj,dt,desiredHeading,desiredSpeed)
            % This function is designed to take a desired heading and
            % velocity and compare them to the kinematic contraints of the
            % agent.
            
            currentVelocity = obj.localState(7:9);
            % CALCULATE THE DESIRED HEADING RATE
            dH = desiredHeading - zeros(3,1);
            Hdot = dH/dt;
            % GET THE ALLOWABLE ANGULAR VELOCITY
            [bounded_Hdot] = obj.boundValue(Hdot,-obj.angularVelocityLimits,obj.angularVelocityLimits);
            % GET THE ACHIEVED HEADING WITH THE 
            achievedHeading = zeros(3,1) + bounded_Hdot*dt;
            
            % CALCULATE THE DESIRED SPEED CHANGE
            currentSpeed = norm(currentVelocity);
            dSpeed = desiredSpeed - currentSpeed; % The proposed change in speed
            acceleration = dSpeed/dt;
            % GET THE ALLOWABLE ACCELERATIONS
            [bounded_a] = obj.boundValue(acceleration,-obj.linearAccelerationLimits(1),obj.linearAccelerationLimits(1)); 
            % GET THE ACHIEVED SPEED
            achievedSpeed = currentSpeed + bounded_a*dt;
            % BOUND THE ABSOLUTE LINEAR VELOCITY
            [achievedSpeed] = obj.boundValue(achievedSpeed,-obj.linearVelocityLimits(1),obj.linearVelocityLimits(1)); 
        end
        
        % ////////////// LOCAL(AGENT) CORE OPERATIONS /////////////////////
        % ALL IN ONE AGENT-INFORMATION UPDATE (FROM ENVIRONMENT)
        function [obj,obstacleSet,agentSet,waypointSet] = getAgentUpdate(obj,observedObjects)
            % This function computes the agents process for updating its
            % knowledge of the environment.
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
            
            for item = 1:length(observedObjects)
                % DEFINE THE SENSED VARIABLES
                objectName = observedObjects(item).name;
                objectID   = observedObjects(item).objectID;               % The objects ID
                objectType = observedObjects(item).type;                   % The objects type (obstacle or waypoint)
                objectSize = observedObjects(item).size;
                objectState = [observedObjects(item).position;observedObjects(item).velocity];
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
            [obj,~] = obj.waypointUpdater(waypointSet);  % Update the target waypoint and heading vector
        end
        % SELECT WAYPOINT
        function [obj,waypointVector] = waypointUpdater(obj,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if ~exist('waypointSet','var') || isempty(waypointSet)
                obj.targetWaypoint = [];       % No waypoints are available
                waypointVector = [1;0;0];      % Default heading
                return
            end
            
            % WAYPOINT SET IS POPULATED -> GET THE PRIORITY VECTOR
            priorityVector = [waypointSet.priority];                       % The vector of priorities
            waypointIndex = linspace(1,length(priorityVector),length(priorityVector)); % Memory of their position in the waypoint set
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
%                     validWaypointIDs = availableWaypointIDs(availableWaypointIDs ~= invalidWaypointIDs);   % Get those visible IDs that are valid
                    validWaypoints = availableWaypointIDs ~= invalidWaypointIDs;   % Get those visible IDs that are valid
                    % IF NO FURTHER VALID IDs
                    if ~any(validWaypoints)
                        obj.targetWaypoint = [];
                        waypointVector = [1;0;0];
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
            waypointAchievedCondition = norm(obj.targetWaypoint.state(1:3)) - (obj.targetWaypoint.size + obj.VIRTUAL.size);
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
                   waypointVector = [1;0;0];
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
            waypointVector = obj.targetWaypoint.state(1:3)/norm(obj.targetWaypoint.state(1:3));  
        end
        % UPDATE OBSTACLE KNOWLEDGE
        function [obj] = updateAgentKnowledge(obj,newEntry)
            % Updates the agents knowledge based on the new information
            % sensed by the agent. 
            % INPUTS:
            % objectID   - The objects ID
            % objectType - The objects sim type (obstacle or waypoint)
            % size       - The objects apparent size
            % newState   - The objects estimated state
            % priority   - The objects priority
            
            % INPUT HANDLING
            if ~isfield(newEntry,'priority')
                newEntry.priority = 0;
            end

            % CHECK CURRENT KNOWLEDGE STATUS
            if isempty(obj.memory)
                obj.memory = vertcat(obj.memory,newEntry); % No knowledge at all
                return
            end
            
            % GET THE LOGICAL INDICES OF ID OCCURANCES
            logicalIDIndex = [obj.memory.objectID] == newEntry.objectID;   % Appearance of object ID in memory
            IDoccurances = sum(logicalIDIndex);
            
            % CHECK FOR THE APPEARANCE OF THE OBJECT IN MEMORY
            if IDoccurances > 1
                error('MEMORY DUPLICATE DETECTED.');
            elseif IDoccurances == 0
                % APPEND ENTRY FOR NEW OBJECT
                obj.memory = vertcat(obj.memory,newEntry); % Append to knowledge to other records
            else
                % UPDATE THE MEMORY STRUCTURE WITH THE NEW ENTRY
                obj.memory(logicalIDIndex) = newEntry;                     % Update existing knowledge of  that agent
            end
        end
        % GET OBJECT MEMORY ITEM FOR OBJECT-ID 
        function [lastEntry] = getLastEntryFromMemory(obj,objectID)
            % This function is designed to return the complete memory item
            % for a given object ID
            % INPUT:
            % obj.memory - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % memoryItem - The agents last memory record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.memory) || ~any([obj.memory.objectID] == objectID)
                lastEntry = [];
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.memory.objectID] == objectID;
            lastEntry = obj.memory(IDindex); 
        end
        % GET STATE FOR OBSTACLE ID 
        function [lastState] = getLastStateFromMemory(obj,objectID)
            % This function is designed to return the last known state
            % vector for given obstacle object ID.
            % INPUT:
            % obj.memory - The agents knowledge matrix 
            % objectID            - The obstacles ID field
            % OUTPUT:
            % lastState - The agents last state record of that ID
                       
            % CASE - ZERO KNOWLEDGE
            if isempty(obj.memory) || ~any([obj.memory.objectID] == objectID)
                lastState = [];
                return
            end
            % EXTRACT THE LAST KNOWN STATE OF THE OBSTACLE 
            IDindex = [obj.memory.objectID] == objectID;
            lastState = obj.memory(IDindex).state; 
        end
        
        % /////////////// OUTPUT FUNCTIONS ////////////////////////////////
        % GET ANOTHER FIGURE FRAME 
        function [obj] = getAnimationFrame(obj,overHandle,TIME,fileName)
            % This function adds the figure to a succession of plots.
            
            figureDimensions = [30 50 800 800];
            set(overHandle, 'Position', figureDimensions); % [x y width height]
            
            % ANNOTATION WITH SIMTIME
            annoString = sprintf('Time: %ss',num2str(TIME.currentTime));
            dim = [0.85 0.05 0.12 0.03];
            annotation('textbox',dim,'String',annoString,'FitBoxToText','off');
            
            % DRAW AND COLLATE GIF FRAME
            im = frame2im(getframe(overHandle));
            [imind,colourMap] = rgb2ind(im,256);
            
            % APPEND THE FRAMES
            slomoFactor = 2;
            if ~obj.giffOn
                imwrite(imind,colourMap,strcat(fileName,'.gif')...
                    ,'gif', 'Loopcount',inf,'DelayTime',TIME.dt*slomoFactor);
                obj.giffOn = 1;
            else
                imwrite(imind,colourMap,strcat(fileName,'.gif')...
                    ,'gif','WriteMode','append','DelayTime',TIME.dt*slomoFactor);
            end
        end
    end
    %% STATIC AGENT TOOLS
    methods (Static)
        % DEFINE THE PITCH AND YAW TO MEET TARGET HEADING
        function [theta,lambda] = getControlInputs(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS
            % V - The current unity velocity vector
            % U - The unit correction vector
            
            % GET THE LINE OF SIGHT ANGLE             
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];                                            % Reject the vertical elements
            % DEFINE LOS ANGLE            
            rotationAxis = cross(Vh,Uh);
            lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            % GET THE ELEVATION ANGLE 
            theta = atan2(U(3),norm(Uh));
        end
        % CALCULATE THE NEW STATE ESTIMATE
        function [newState] = linearStateEstimation(priorState,newPosition,dt)
           % This function takes the previous known state of the obstacle
           % and estimates its new state.
           % INPUTS:
           % priorState  - [Sx;Sy;Sz;Vx;Vy;Vz];
           % newPosition - The new relative position
           % OUTPUT:
           % newState - The new obstacle state
 
           % GET THE PRIOR VALUES
           position0 = priorState(1:3);
           velocity0 = priorState(4:6);
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
        
        % GET ESTIMATE OF RELATIVE POSITION & VELOCITY
        function [cartesianPosition] = cartesianFromSpherical(range,azimuth,elevation)
           % This function calculates the relative position from a spherical 
           % coordinate system that assumes angles are measured from the 
           % agents heading vector Xref.
           % INPUTS:
           % range     - The objects radial seperation (m)
           % azimuth   - The objects relative azimuth angle (rad)
           % elevation - The objects relative elevation (rad)
           % OUTPUTS:
           % cartesianPosition  - The new 3D position in local coordinates (m)

           % DEFINE THE POSITION AS A VECTOR INTERVAL
           cartesianPosition = [cos(azimuth)*cos(elevation);...
                                sin(azimuth)*cos(elevation);...
                                             sin(elevation)]*range;
        end
        % CONVERT CARTESIAN TO SPHERICAL
        function [range,azimuth,elevation] = sphericalFromCartesian(positionVector)
            range = norm(positionVector);                                            % The range
            azimuth = atan2(positionVector(2),positionVector(1));                    % The angle made in the azimuth (bearing)
            elevation = atan2(positionVector(3),sqrt(positionVector(1).^2 + positionVector(2).^2));  % The elevation (vertical bearing)
        end
        % VALIDATE THE OBSTACLE
        function [tau] = validateCollision(S_b,c_b)
            % CONFIRM WE KNOW THE HEADING OF THE OBSTACLE
            if any(isnan(c_b))
                tau = -inf;
                return
            end                 
            % Define the time to closest approach (+ve converging)
            tau = -(dot(S_b,c_b)/dot(c_b,c_b));
        end
    end  
end

