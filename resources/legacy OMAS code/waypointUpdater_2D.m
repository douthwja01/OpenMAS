% SELECT WAYPOINT (FOR 2D STATES)
function [obj,waypointVector] = waypointUpdater_2D(obj,waypointSet)
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
waypointAchievedCondition = norm(obj.targetWaypoint.position) - (obj.targetWaypoint.radius + obj.VIRTUAL.radius);
% IF THE CRITERIA IS MET AND THE ID IS NOT LOGGED
if waypointAchievedCondition <= 0 && ~any(ismember(obj.achievedWaypoints,obj.targetWaypoint.objectID))
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
waypointVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
end