%% AGENT TOOLS (agent_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class aims to abstract the simulation/technical elements of agents
% definitions.

classdef agent_tools < objectDefinition
    
    properties
        % OBSTACLE OBJECT MATRIX (record of sightings)
        MEMORY;                         % Record of relative scene
        maxSamples = 1;                 % The maximum number of states retained
    end
            
    % /////////////////////// SIMULATION INTERFACES ///////////////////////
    methods
        % UPDATE OBJECTIVE
        function [obj,waypointVector] = UpdateTargetWaypoint(obj,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if ~exist('waypointSet','var') || isempty(waypointSet)
                obj.targetWaypoint = [];                % No waypoints are available
                waypointVector = [];                    % Default heading
                obj = obj.SetIdleStatus(logical(true)); % Set idle
                return
            end
            
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
                        obj = obj.SetIdleStatus(logical(true));            % Set the idle status
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
            waypointGetCondition = obj.GetTargetCondition();
            
            % IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
            if waypointGetCondition && ~any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID))
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
            [waypointVector] = obj.GetTargetHeading();
        end
        % GET TARGET HEADING VECTOR OF VISIBLE OBJECT
        function [headingVector] = GetTargetHeading(obj,targetObject)
            % This function calculates the heading vector to the current
            % obj.targetWaypoint, or to the provided object with a position
            % field.
            
            % If a target is provided
            if nargin > 1
                targetPosition = obj.GetLastMeasurement(targetObject.objectID,'position');
            elseif ~isempty(obj.targetWaypoint)
                targetPosition = obj.targetWaypoint.position(:,obj.targetWaypoint.sampleNum);
            else
                % Default heading
                if obj.VIRTUAL.is3D
                    targetPosition = [1;0;0];
                else
                    targetPosition = [1;0];
                end
            end
            % Get the vector heading of the target
            headingVector = targetPosition/norm(targetPosition);
        end
        % GET TARGET-GET CONDITION
        function [targetLogical] = GetTargetCondition(obj)
            % Get the current measurements
            currentPosition = obj.targetWaypoint.position(:,obj.targetWaypoint.sampleNum);
            currentRadius   = obj.targetWaypoint.radius(obj.targetWaypoint.sampleNum);
            % Check the current achieve-tolerances on a given target
            targetLogical = 0 > norm(currentPosition) - (currentRadius + obj.VIRTUAL.radius);
        end
    end
    % MEMORY OPERATIONS
    methods
        % //////////////////// MEMORY MANIPULATION ////////////////////////
        % SORT MEMORY - DEFINED FIELD
        function [obj] = SortMemoryByField(obj,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            assert(isfield(obj.MEMORY,field),'Field must belong to the memory structure.');
            
            % Reorder the memory structure based on fieldname
            [~,ind] = sort([obj.MEMORY.(field)],2,'descend');         % Ordered indices of the object IDs
            % Sort the memory structure 
            obj.MEMORY = obj.MEMORY(ind);
        end
        % GET OBJECT PRIORITY
        function [priority_j] = GetObjectPriority(obj,objectID)
            logicalIDIndex = [obj.MEMORY(:).objectID] == objectID;
        	priority_j = 1/norm(obj.MEMORY(logicalIDIndex).position(:,obj.MEMORY(logicalIDIndex).sampleNum));        
        end
        % Set the maximum sample number
        function [obj] = SetBufferSize(obj,horizon)
            assert(numel(horizon) == 1 && isnumeric(horizon),'Horizon must be a scalar number of steps.');
            % define the size of the buffer
            obj.maxSamples = horizon;
            % Initialise memory structure (with 3D varient)
            obj.MEMORY = obj.GetMemoryStructure(obj.maxSamples);
        end
            
        % ///////////////////////// MEMORY I/O ////////////////////////////
        % MEMORY - SET LAST MEASUREMENT BY FIELD
        function [obj]  = SetLastMeasurementByObjectID(obj,objectID,field,data)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            
            % Check occurances
            if ~obj.IsInMemory(objectID)
                return      % No data to return
            end
            
            % Check the field exists in memory
            assert(isfield(obj.MEMORY,field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            
            % Get the associated buffer data
            memData = obj.MEMORY(memoryLogicals).(field);
            
            
            
            % Get the latest samples data from buffer  
            if isa(memData,'circularBuffer')
                memData(:,obj.MEMORY(memoryLogicals).sampleNum) = data;
            else
                memData(:,1) = data;
            end
            % Reapply the data to the memory field [TO-DO] Rewrite
            obj.MEMORY(memoryLogicals).(field) = memData;
        end
        % MEMORY - UPDATE FROM OBSERVATIONS
        function [obj]  = UpdateMemoryFromObject(obj,t_sample,observedObject)
            
            % This program is designed to parse the sensor data recieved into memory
            % for access in main agent code.
            
            % Input sanity check #1
            if numel(observedObject) == 0
                return
            end
                       
            % [TO-DO] REMOVE MEMORY ENTRIES NOT VISIBLE ///////////////////
            % Get the indices of objects still visible
            % visibleIDs = ismember([obj.MEMORY.objectID],[observedObject.objectID]);
            % Remove entries no longer visible
            % obj.MEMORY = obj.MEMORY(visibleIDs);                              % Removes all entries not visible 
            
            % ////////// UPDATE MEMORY WITH NEW ENTRIES AND EXISTING ENTRIES //////////
            % Get the logical indices of ID occurances
            logicalIDIndex = [obj.MEMORY.objectID] == observedObject.objectID;  % Appearance of object ID in memory
            IDoccurances   = sum(logicalIDIndex);
            % Determine memory behaviour
            switch IDoccurances
                case 0 % ///////////////// IS NOT IN MEMORY YET ///////////////////
                    % Override the template if its the first reading
                    if numel(obj.MEMORY) == 1 && obj.MEMORY(1).objectID == 0
                        logicalIDIndex = 1;                                     % The memory address to be overwritten
                    else
                        % Get new memory structure with associated fields
                        newEntry = obj.GetMemoryStructure(obj.maxSamples);
                        obj.MEMORY = [obj.MEMORY;newEntry];                     % Append the new structure
                        logicalIDIndex = numel(obj.MEMORY);                     % Address is the end address
                    end
                    
                    % Record current time and append to sample 
                    obj.MEMORY(logicalIDIndex).time(1) = t_sample;
                    
                    % UPDATE MEMORY FIELDS FROM OBSERVATIONS
                    memFields = fieldnames(obj.MEMORY(logicalIDIndex));         % Update fields by dynamic association
                    for entry = 1:numel(memFields)
                        if ~isfield(observedObject,memFields{entry})            % Check the memory field is observable
                            continue                                            % Field cannot be updated from observations
                        end
                        % Data is available
                        if isa(obj.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            % Update circular-buffers
                            data = obj.MEMORY(logicalIDIndex).(memFields{entry});
                            data(:,obj.MEMORY(logicalIDIndex).sampleNum) = observedObject.(memFields{entry});
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                        else
                            % Direct value override
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = observedObject.(memFields{entry});
                        end
                        
                    end
                case 1 % ////////// IS IN MEMORY ALREADY AND SINGUAR //////////////
                    % Indicate new sample
                    obj.MEMORY(logicalIDIndex).sampleNum = obj.MEMORY(logicalIDIndex).sampleNum + 1;
                    % Get the new sample number
                    n = obj.MEMORY(logicalIDIndex).sampleNum; 
                    
                    % Record current time and append to sample 
                    obj.MEMORY(logicalIDIndex).time(:,n) = t_sample;

                    % Update dynamic fields by association
                    memFields = fieldnames(obj.MEMORY(logicalIDIndex));
                    for entry = 1:numel(memFields)
                        if ~isfield(observedObject,memFields{entry})            % Check the memory field is observable
                            continue                                            % Field cannot be updated from observations
                        end
                        % Isolate circular-buffers
                        if isa(obj.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            data = obj.MEMORY(logicalIDIndex).(memFields{entry});
                            data(:,n) = observedObject.(memFields{entry});
                            obj.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                        end
                    end
                otherwise
                    error('[ERROR] Agent memory structure is distorted.');
            end
            % /////////// DEDUCTIONS FROM THE MEASUREMENTS ////////////
            % Call external priority calculator (based on range, position
            % speed.. etc)
            obj.MEMORY(logicalIDIndex).priority = obj.GetObjectPriority(observedObject.objectID);
        end
        % MEMORY - WHOLE STRUCTURE BY ID
        function [data] = GetMemoryStructByID(obj,objectID)
            % Check occurances
            data = [];
            if ~obj.IsInMemory(obj,objectID)
                return
            end
            % Get the location of data in the memory structure
            data = obj.MEMORY([obj.MEMORY(:).objectID] == objectID);                        % The memory structure
        end
        % MEMORY - TRAJECTORY BY FIELD
        function [data] = GetTrajectoryByObjectID(obj,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            assert(ischar(field),'Field must be a defined as a string label.');
            % Check occurances
            data = [];
            if ~obj.IsInMemory(objectID)
                return      % No data to return
            end
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            data = obj.MEMORY(memoryLogicals).(field);
            % Unpack the buffer along the [:,n] dimension
            data = data.unpack(2,obj.MEMORY(memoryLogicals).sampleNum);
        end
        % MEMORY - LAST MEASUREMENT BY FIELD
        function [data] = GetLastMeasurementByObjectID(obj,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            
            % Check occurances
            data = [];
            if ~obj.IsInMemory(objectID)
                return      % No data to return
            end
            
            % Check the field exists in memory
            assert(isfield(obj.MEMORY,field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            data = obj.MEMORY(memoryLogicals).(field);
            % Get the latest samples data from buffer  
            if isa(data,'circularBuffer')
                data = data(:,obj.MEMORY(memoryLogicals).sampleNum);
            else
                data = data(:,1);
            end
        end
        % MEMORY - CHECK ID IS PRESENT
        function [flag] = IsInMemory(obj,objectID)
            % Get the location of data in the memory structure
            memoryLogicals = [obj.MEMORY(:).objectID] == objectID;
            flag = any(memoryLogicals); % Check for occurances
        end
    end
    methods (Static)
        % //////////////////// MEMORY INITIALISATION //////////////////////
        % FROM MEMORY - LAST MEASUREMENT FROM STRUCTURE
        function [latestData] = GetLastMeasurementFromStruct(memStruct,field)
            % Input sanity check
            assert(ischar(field),'Field must be a defined as a string label.');
            
            % Get the associated buffer data
            latestData = memStruct.(field);
            % Get the latest samples data from buffer  
            if isa(latestData,'circularBuffer')
                latestData = latestData(:,memStruct.sampleNum);
            else
                latestData = latestData(:,1);
            end
        end
        % GET EMPTY MEMORY STRUCTURE (3D trajectories)
        function [memStruct]  = GetMemoryStructure(horizonSteps)
            % This function contains a basic agent-memory structure. This
            % is used to retain information on observed objects and maintain 
            % a regular structure.
            
            % Input sanity check
            if nargin < 1
                horizonSteps = 10; % Duration retained in memory
            end
            
            % The fields of memory structure define the fields of the
            % simulation 'observation' structure that are retained.
            
            % Create empty memory structure
            memStruct = struct(...
                'name','',...
                'objectID',uint8(0),...
                'type',OMAS_objectType.misc,...
                'sampleNum',uint8(1),...
                'time',circularBuffer(NaN(1,horizonSteps)),...
                'position',circularBuffer(NaN(3,horizonSteps)),...
                'velocity',circularBuffer(NaN(3,horizonSteps)),...
                'radius',circularBuffer(NaN(1,horizonSteps)),...
                'range',circularBuffer(NaN(1,horizonSteps)),...
                'heading',circularBuffer(NaN(1,horizonSteps)),...
                'elevation',circularBuffer(NaN(1,horizonSteps)),...
                'width',circularBuffer(NaN(1,horizonSteps)),...
                'geometry',struct('vertices',[],'faces',[],'normals',[],'centroid',[]),...
                'priority',[]);
        end
    end
end