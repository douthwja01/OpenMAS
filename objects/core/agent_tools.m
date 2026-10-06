%% AGENT TOOLS (agent_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class aims to abstract the simulation/technical elements of agents
% definitions.

classdef agent_tools < objectDefinition
    %% GENERIC PROPERTIES
    properties
        % OBSTACLE OBJECT MATRIX (record of sightings)
        MEMORY;           % Record of relative scene
        maxSamples = 1;   % The maximum number of states retained
    end
    %% SIMULATION INTERFACES
    methods
        % UPDATE OBJECTIVE
        function [this,waypointVector] = UpdateTargetWaypoint(this,waypointSet)
            % This function returns the heading for the waypoint with the
            % highest listed priority. The agent retains a matrix of
            % objectID's for the waypoints it has achieved. This is used to
            % select the next priority waypoint.
            
            % INPUT HANDLING
            if isempty(waypointSet)
                this.targetWaypoint = [];           % No waypoints are available
                waypointVector = [];                % Default heading
                this.SetGLOBAL('idleStatus',true); 	% Set idle
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
            if isempty(this.targetWaypoint)
                % NO WAYPOINTS HAVE BEEN SELECTED
                if isempty(this.achievedWaypoints)
                    [~,maxPriorityIndex] = max(waypointMatrix(3,:));       % Index of maximum priority value in waypoint matrix
                    waypointSetIndex = waypointMatrix(1,maxPriorityIndex); % Gets the associated waypoint set index of the highest priority
                    this.targetWaypoint = waypointSet(waypointSetIndex);    % Update with the highest priority waypoint
                else
                    % SELECT HIGHEST PRIORITY, BUT NOT ACHIEVED
                    invalidWaypointIDs = this.achievedWaypoints;            % Get the achieved waypoint IDs
                    availableWaypointIDs = waypointMatrix(2,:);            % Get those visable IDs
                    validWaypoints = availableWaypointIDs ~= invalidWaypointIDs;   % Get those visible IDs that are valid
                    % IF NO FURTHER VALID IDs
                    if ~any(validWaypoints)
                        this.targetWaypoint = [];
                        waypointVector = [];
                        this.SetGLOBAL('IdleStatus',true);                  % Set the idle status
                        return
                    end
                    
                    % SELECT REMAINING WAYPOINTS
                    validWaypoints = waypointMatrix(:,validWaypoints); % Get waypoinSet location, IDs and priority of valid waypoints
                    % SELECT WAYPOINT OF NEXT HIGHEST PRIORITY
                    [~,maxValidIndex ] = max(validWaypoints(3,:));       % Get the index of max priority waypoint
                    waypointSetIndex = validWaypoints(1,maxValidIndex);  % Get the location of the target waypoint in the waypoint set
                    this.targetWaypoint = waypointSet(waypointSetIndex);  % Select waypoint object
                end
            end
            
            % EVALUATE CURRENT TARGET WAYPOINT ////////////////////////////
            % CHECK THE TOLERANCES ON THE CURRENT TARGET WAYPOINT
            waypointGetCondition = this.GetTargetCondition();
            
            % IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
            if waypointGetCondition && ~any(ismember(this.achievedWaypoints,this.targetWaypoint.objectID))
                this.achievedWaypoints = horzcat(this.achievedWaypoints,this.targetWaypoint.objectID); % Add achieved objectID to set
            end
            
            % CHECK IF CURRENT TARGET IS VALID ////////////////////////////
            % ASSESS THE ACHIEVED WAYPOINT SET AGAINST THE CURRENT TARGET
            invalidTargetCondition = any(ismember(this.achievedWaypoints,this.targetWaypoint.objectID)); % Is the current target ID in the achieved set
            if invalidTargetCondition
                % TARGET WAYPOINT HAS BEEN ACHIEVED, REALLOCATE ONE
                selectionVector = priorityVector < this.targetWaypoint.priority;   % Only those with lower priority
                if ~any(selectionVector)
                    this.targetWaypoint = [];
                    waypointVector = [];
                    return
                else
                    reducedMatrix = waypointMatrix(:,selectionVector);            % Build vector of priorities less than
                    [~,reducedIndex] = max(reducedMatrix(3,:));                   % Finds the max priority index where less than current
                    waypointIndex = reducedMatrix(1,reducedIndex);                % Get the waypoint index from the reduced priority matrix
                    this.targetWaypoint = waypointSet(waypointIndex);              % Select the next highest priority waypoint
                end
            else
                currentWaypointID = this.targetWaypoint.objectID;                 % Get the ID of the current target waypoint
                selector = (waypointMatrix(2,:) == currentWaypointID);           % Find where the ID appears in the waypoint matrix
                waypointIndex = waypointMatrix(1,selector);                      % Get the waypoint Set index from the waypoint matrix
                this.targetWaypoint = waypointSet(waypointIndex);                 % Update with the highest priority waypoint
            end
            % THE CORRECT TARGET WAYPOINT IS NOW DEFINED, GENERATE HEADING
            [waypointVector] = this.GetTargetHeading();
        end
        % GET TARGET HEADING VECTOR OF VISIBLE OBJECT
        function [headingVector] = GetTargetHeading(this,targetObject)
            % This function calculates the heading vector to the current
            % this.targetWaypoint, or to the provided object with a position
            % field.
            
            % If a target is provided
            if nargin > 1
                targetPosition = this.GetLastMeasurementByID(targetObject.objectID,'position');
            elseif ~isempty(this.targetWaypoint)
                targetPosition = this.targetWaypoint.position(:,this.targetWaypoint.sampleNum);
            else
                % Default heading
                if this.Is3D
                    targetPosition = [1;0;0];
                else
                    targetPosition = [1;0];
                end
            end
            % Get the vector heading of the target
            headingVector = targetPosition/norm(targetPosition);
        end
        % GET TARGET-GET CONDITION
        function [targetLogical] = GetTargetCondition(this)
            % Get the current measurements
            waypointPosition = this.targetWaypoint.position(:,this.targetWaypoint.sampleNum);
            waypointRadius   = this.targetWaypoint.radius(this.targetWaypoint.sampleNum);
            radialConstraint = this.GetGLOBAL('radius');
            % Check the current achieve-tolerances on a given target
            targetLogical = 0 > (norm(waypointPosition) - (waypointRadius + radialConstraint));
        end
    end
    %% UPDATE METHODS
    methods
        % GLOBAL UPDATE - EULER 6DOF(3DOF) TO NED STATE VECTOR
        function [this] = GlobalUpdate_NED(this,dt,eulerState)
            
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
                eulerIndices = 4:6;
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
                eulerIndices = 3;
            else
                error('State notation not recognised');
            end
            
            % Previous state 'k' 
            position_k   = this.GetGLOBAL('position');
            velocity_k   = this.GetGLOBAL('velocity');
            quaternion_k = this.GetGLOBAL('quaternion');
            state_k      = this.GetGLOBAL('priorState');
            
            % EQUIVALENT RATES
            velocity_k_plus   = (eulerState(positionIndices) - state_k(positionIndices))/dt;
            eulerRates_k_plus = (eulerState(eulerIndices) - state_k(eulerIndices))/dt;
            
            % NED ROTATION CONVENTION
            % 3D - Rotation order (Z Y X)
            % 2D - Rotation order (Z)
            
            % Formulate the conversion
            globalAxisRates = eulerRates_k_plus;
            if numel(globalAxisRates) ~= 3
                globalAxisRates = zeros(3,1) + [-1;0;0]*eulerRates_k_plus;
                velocity_k_plus = [velocity_k_plus;0];
            else
                % GET CONSTANT CONVERSION FROM GLOBAL NED TO GLOBAL ENU
                R_NED_ENU = [1 0 0 ; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
                % Formulate the conversion
                globalAxisRates = W*eulerRates_k_plus;
            end
            
            % INTEGRATE THE GLOBAL QUATERNION
            [quaternion_k_plus] = OMAS_geometry.integrateQuaternion(quaternion_k,globalAxisRates,dt);
            % NEW ROTATION MATRIX (G>B)
            R_k_plus = OMAS_geometry.quaternionToRotationMatrix(quaternion_k_plus);
            
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = position_k + dt*velocity_k;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.GlobalUpdate_direct(...
                globalPosition_k_plus,... % Global position at k plius
                globalVelocity_k_plus,... % Global velocity at k plus
                quaternion_k_plus,...     % Quaternion at k plus
                eulerState);              % The new state for reference
        end
        % GLOBAL UPDATE - EULER 6DOF(3DOF) TO NED STATE VECTOR ( LOCAL AXIS REMAINS STATIC )
        function [this] = GlobalUpdate_NED_fixed(this,dt,eulerState)
            
            % NOTATION INDICIES
            if numel(eulerState) == 6
                positionIndices = 1:3;
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
            else
                error('State notation not recognised');
            end
            
            % Previous state 'k' 
            position_k   = this.GetGLOBAL('position');
            velocity_k   = this.GetGLOBAL('velocity');
            quaternion_k = this.GetGLOBAL('quaternion');
            state_k      = this.GetGLOBAL('priorState');
            
            % EQUIVALENT RATES
%             velocity_k_plus   = (eulerState(positionIndices) - state_k(positionIndices))/dt;
            
            % NED ROTATION CONVENTION
            % 3D - Rotation order (Z Y X)
            % 2D - Rotation order (Z)
            % ROTATION RATES ABOUT THE GLOBAL AXES
            
            % /////////////// DEFINE PREVIOUS PARAMETERS //////////////////
            velocity_k_plus = eulerState(velocityIndices,1);    % The velocity in the local NED frame
            
            % Formulate the conversion
            if numel(velocity_k_plus) ~= 3
                velocity_k_plus = [velocity_k_plus;0];
            end
            
            % GET LOCAL NED TO GLOBAL NED ROTATION MATRIX
            R_k = OMAS_geometry.quaternionToRotationMatrix(quaternion_k);
            % GET CONSTANT CONVERSION FROM GLOBAL NED TO GLOBAL ENU
            R_NED_ENU = [1 0 0 ; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
            
            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            velocity_k_plus = R_NED_ENU*(R_k*velocity_k_plus);
            position_k_plus = position_k + dt*velocity_k;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [this] = this.GlobalUpdate_direct(...
                position_k_plus,...     % Global position at k plius
                velocity_k_plus,...     % Global velocity at k plus
                quaternion_k_plus,...   % Quaternion at k plus
                eulerState);            % The new state for reference
        end
        % GLOBAL UPDATE - DIRECT (AGENT OVERRIDE)
        function [this] = GlobalUpdate_direct(this,p,v,q)
            % Under this notation, the state vector already contains the
            % global parameters of the object.
            % INPUTS:
            % globalPosition - 3D global cartesian position given in the ENU coordinate frame.
            % globalVelocity - 3D global cartesian velocity given in the ENU coordinate frame.
            % quaternion     - The new quaternion pose of body in the global frame.
            % R              - The rotation of the body
            % this.localState - The previous localState (independant of convention)
            % eulerState     - The new state as reported by the agent
            
            % Input sanity check
            assert(IsColumn(p,3),'Global position must be a 3D column vector [3x1].');
            assert(IsColumn(v,3),'Global velocity must be a 3D column vector [3x1].');
            assert(IsColumn(q,4),'Global pose must be a 4D quaternion vector [4x1].');
            assert(size(this.localState,2) == 1,'The length of the objects state update must match the its local state.');
            
            % Check if the object idle condition is made
            if this.GetGLOBAL('type') == OMAS_objectType.agent
                % Evaluate way-point
                if isempty(this.targetWaypoint) && ~isempty(this.achievedWaypoints)
                    this = this.SetGLOBAL('idleStatus',true);
                end
            end
            
            % ///////////////// REASSIGN K+1 PARAMETERS //////////////////////////
            % Convert the quaternion to the equivalent rotation matrix
            R = OMAS_geometry.quaternionToRotationMatrix(q);
            % Assign the global parameters
            this.SetGLOBAL('position',p);                 	% Reassign the global position
            this.SetGLOBAL('velocity',v);                 	% Reassign the global velocity
            this.SetGLOBAL('quaternion',q);                	% Reassign the quaternion
            this.SetGLOBAL('R',R);                          % New rotation matrix
            this.SetGLOBAL('priorState',this.localState);  	% Record the previous state
        end
    end
    % Only the objectDefinition class has access
    methods (Static)
        % Create agent VIRTUAL structure (Override the objectDefinition)
        function [GLOBAL]  = CreateGLOBAL()
            % This is function that assembles the template desciption of an
            % agent type object.
            
            % Define the VIRTUAL structure
            GLOBAL = struct();
            GLOBAL.type = OMAS_objectType.agent;
            GLOBAL.hitBoxType = OMAS_hitBoxType.spherical;
            GLOBAL.detectionRadius = inf;
            GLOBAL.radius = 0.5;                % Diameter of 1m
            GLOBAL.colour = rand(1,3,'single');	% Colour (for plotting)
            GLOBAL.symbol = 'diamond';          % Representative symbol
            GLOBAL.position = [0;0;0];          % Global Cartesian position
            GLOBAL.velocity = [0;0;0];          % Global Cartesian velocity
            GLOBAL.quaternion = [1;0;0;0];      % Global quaternion pose
            GLOBAL.idleStatus = false;          % Object idle logical
            GLOBAL.is3D = true;                 % Object operates in 3D logical
            GLOBAL.priorState = [];
        end
    end
    %% MEMORY MANIPULATION
    methods
        % MEMORY - Sort by defined field
        function [this] = UpdateMemoryOrderByField(this,field)
            % INPUTS:
            % type - Sort option; memory field label.
            
            % Input sanity check
            assert(ischar(field),'Memory sort method must be a string.');
            assert(isfield(this.MEMORY,field),'Field must belong to the memory structure.');
            
            % Reorder the memory structure based on fieldname
            [~,ind] = sort([this.MEMORY.(field)],2,'descend');         % Ordered indices of the object IDs
            % Sort the memory structure
            this.MEMORY = this.MEMORY(ind);
        end
        % MEMORY - Update memory from an observation structure
        function [this] = UpdateMemoryFromObject(this,t_sample,observedObject)
            
            % This program is designed to parse the sensor data recieved into memory
            % for access in main agent code.
            
            % Input sanity check #1
            if numel(observedObject) == 0
                return
            end
            
            % [TO-DO] REMOVE MEMORY ENTRIES NOT VISIBLE ///////////////////
            % Get the indices of objects still visible
            % visibleIDs = ismember([this.MEMORY.objectID],[observedObject.objectID]);
            % Remove entries no longer visible                      
            
            % ////////// UPDATE MEMORY WITH NEW ENTRIES AND EXISTING ENTRIES //////////
            % Get the logical indices of ID occurances
            logicalIDIndex = [this.MEMORY.objectID] == observedObject.objectID;  % Appearance of object ID in memory
            IDoccurances   = sum(logicalIDIndex);
            % Determine memory behaviour
            switch IDoccurances
                case 0 % ///////// IS NOT IN MEMORY YET ///////////////////
                    % Override the template if its the first reading
                    if numel(this.MEMORY) == 1 && this.MEMORY(1).objectID == 0
                        logicalIDIndex = 1;                                     % The memory address to be overwritten
                    else
                        % Get new memory structure with associated fields
                        newEntry = this.CreateMEMORY(this.maxSamples);
                        this.MEMORY = [this.MEMORY;newEntry];                     % Append the new structure
                        logicalIDIndex = numel(this.MEMORY);                     % Address is the end address
                    end
                    
                    % Record current time and append to sample
                    this.MEMORY(logicalIDIndex).time(1) = t_sample;
                    
                    % UPDATE MEMORY FIELDS FROM OBSERVATIONS
                    memFields = fieldnames(this.MEMORY(logicalIDIndex));         % Update fields by dynamic association
                    for entry = 1:numel(memFields)
                        if ~isfield(observedObject,memFields{entry})            % Check the memory field is observable
                            continue                                            % Field cannot be updated from observations
                        end
                        % Data is available
                        if ~isa(this.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            % Direct value override
                            data = observedObject.(memFields{entry});
                        else
                            % Update circular-buffers
                            data = this.MEMORY(logicalIDIndex).(memFields{entry});
                            % Attempt to merge data
                            try
                                data(:,this.MEMORY(logicalIDIndex).sampleNum) = observedObject.(memFields{entry});
                            catch memoryUpdateError
                                %warning('Error inserting new data entry, are they the same dimension?')
                                rethrow(memoryUpdateError);
                            end
                        end
                        % Attempt to apply
                        this.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                    end
                case 1 % ////////// IS IN MEMORY ALREADY AND SINGUAR //////////////
                    % Indicate new sample
                    this.MEMORY(logicalIDIndex).sampleNum = this.MEMORY(logicalIDIndex).sampleNum + 1;
                    % Get the new sample number
                    n = this.MEMORY(logicalIDIndex).sampleNum;
                    
                    % Record current time and append to sample
                    this.MEMORY(logicalIDIndex).time(:,n) = t_sample;
                    
                    % Update dynamic fields by association
                    memFields = fieldnames(this.MEMORY(logicalIDIndex));
                    for entry = 1:numel(memFields)
                        if ~isfield(observedObject,memFields{entry})            % Check the memory field is observable
                            continue                                            % Field cannot be updated from observations
                        end
                        %                         % Isolate circular-buffers
                        %                         if isa(this.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                        %                             data = this.MEMORY(logicalIDIndex).(memFields{entry});
                        %                             data(:,n) = observedObject.(memFields{entry});
                        %                             this.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                        %                         else
                        %
                        %                         end
                        % Data is circular buffer.
                        if isa(this.MEMORY(logicalIDIndex).(memFields{entry}),'circularBuffer')
                            % Update circular-buffers
                            data = this.MEMORY(logicalIDIndex).(memFields{entry});
                            % Attempt to merge data
                            data(:,n) = observedObject.(memFields{entry});
                        else
                            % Direct value override
                            data = observedObject.(memFields{entry});
                        end
                        % Attempt to apply
                        this.MEMORY(logicalIDIndex).(memFields{entry}) = data;
                    end
                otherwise
                    error('[ERROR] Agent memory structure is distorted.');
            end
            % /////////// DEDUCTIONS FROM THE MEASUREMENTS ////////////
            % Call external priority calculator (based on range, position
            % speed.. etc)
            this.MEMORY(logicalIDIndex).priority = this.GetObjectPriority(observedObject.objectID);
        end
        % MEMORY - LAST MEASUREMENT BY FIELD
        function [data] = GetLastMeasurementByID(this,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            
            % Check occurances
            data = [];
            if ~this.IsInMemory(objectID)
                return      % No data to return
            end
            
            % Check the field exists in memory
            assert(isfield(this.MEMORY,field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [this.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            data = this.MEMORY(memoryLogicals).(field);
            % Get the latest samples data from buffer
            if isa(data,'circularBuffer')
                data = data(:,this.MEMORY(memoryLogicals).sampleNum);
            else
                data = data(:,1);
            end
        end
        % MEMORY - SET LAST MEASUREMENT BY FIELD
        function [this] = SetLastMeasurementByID(this,objectID,field,data)
            % This function sets the last value at the given memory field
            % to current value.
            
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            
            % Check occurances
            if ~this.IsInMemory(objectID)
                return      % No data to return
            end
            
            % Check the field exists in memory
            assert(isfield(this.MEMORY,field),'Field must be a defined as a string label.');
            
            % Get the location of data in the memory structure
            memoryLogicals = [this.MEMORY(:).objectID] == objectID;
            
            % Get the associated buffer data
            memData = this.MEMORY(memoryLogicals).(field);
            
            % Get the latest samples data from buffer
            if isa(memData,'circularBuffer')
                memData(:,this.MEMORY(memoryLogicals).sampleNum) = data;
            else
                memData(:,1) = data;
            end
            % Reapply the data to the memory field [TO-DO] Rewrite
            this.MEMORY(memoryLogicals).(field) = memData;
        end
        % MEMORY - Set the maximum sample number
        function [this] = SetBufferSize(this,horizon)
            assert(numel(horizon) == 1 && isnumeric(horizon),'Horizon must be a scalar number of steps.');
            % define the size of the buffer
            this.maxSamples = horizon;
            % Initialise memory structure (with 3D varient)
            this.MEMORY = this.CreateMEMORY(this.maxSamples);
        end
        % MEMORY - Get the object priority
        function [p_j]  = GetObjectPriority(this,objectID)
            logicalIDIndex = [this.MEMORY(:).objectID] == objectID;
            p_j = 1/norm(this.MEMORY(logicalIDIndex).position(:,this.MEMORY(logicalIDIndex).sampleNum));
        end
        % MEMORY - WHOLE STRUCTURE BY ID
        function [data] = GetMemoryStructByID(this,objectID)
            % Check occurances
            data = [];
            if ~this.IsInMemory(this,objectID)
                return
            end
            % Get the location of data in the memory structure
            data = this.MEMORY([this.MEMORY(:).objectID] == objectID);                        % The memory structure
        end
        % MEMORY - TRAJECTORY BY FIELD
        function [data] = GetTrajectoryByID(this,objectID,field)
            % Input sanity check
            assert(isnumeric(objectID),'ObjectID must be a numeric ID.');
            assert(ischar(field),'Field must be a defined as a string label.');
            % Check occurances
            data = [];
            if ~this.IsInMemory(objectID)
                return      % No data to return
            end
            % Get the location of data in the memory structure
            memoryLogicals = [this.MEMORY(:).objectID] == objectID;
            % Get the associated buffer data
            data = this.MEMORY(memoryLogicals).(field);
            % Unpack the buffer along the [:,n] dimension
            data = data.unpack(2,this.MEMORY(memoryLogicals).sampleNum);
        end
        % MEMORY - CHECK ID IS PRESENT
        function [flag] = IsInMemory(this,objectID)
            % Get the location of data in the memory structure
            memoryLogicals = [this.MEMORY(:).objectID] == objectID;
            flag = any(memoryLogicals); % Check for occurances
        end
        % GET EMPTY MEMORY STRUCTURE (2D & 3D trajectories)
        function [MEMORY] = CreateMEMORY(this,horizonSteps)
            % This function contains a basic agent-memory structure. This
            % is used to retain information on observed objects and maintain
            % a regular structure.
            
            % Input sanity check
            if nargin < 1
                horizonSteps = 1;   % Duration retained in memory
            end
            
            % Handle interval memory types (2D & 3D)
            if this.Is3D()
                dim = 3;
            else
                dim = 2;
            end
            
            % The fields of memory structure define the fields of the
            % simulation 'observation' structure that are retained.
            
            % ///// Create empty memory structure /////
            % House-keeping
            MEMORY = struct();
            MEMORY.name = '';
            MEMORY.objectID = uint8(0);
            MEMORY.type = OMAS_objectType.misc;
            MEMORY.sampleNum = uint8(1);
            % Cartesian measurements
            MEMORY.time         = circularBuffer(NaN(1,horizonSteps));
            MEMORY.position     = circularBuffer(NaN(dim,horizonSteps));
            MEMORY.velocity     = circularBuffer(NaN(dim,horizonSteps));
            MEMORY.radius       = circularBuffer(NaN(1,horizonSteps));
            % Spherical measurements
            MEMORY.range        = circularBuffer(NaN(1,horizonSteps));
            MEMORY.heading      = circularBuffer(NaN(1,horizonSteps));
            MEMORY.elevation	= circularBuffer(NaN(1,horizonSteps));
            MEMORY.width    	= circularBuffer(NaN(1,horizonSteps));
            % Additionals
            MEMORY.geometry = struct('vertices',[],'faces',[],'normals',[],'centroid',[]);
            MEMORY.priority = [];
        end
    end
    methods (Static)
        % FROM MEMORY - LAST MEASUREMENT FROM STRUCTURE
        function [latestData] = GetLastMeasurementFromStruct(memStruct,field)
            % Input sanity check
            assert(isstruct(memStruct),'Expecting a valid memory structure.');
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
        
    end
end