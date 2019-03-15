%% OPENMAS TRAJECTORY DATA EXTRACTION (simulation_getTrajectory.m) %%%%%%%%
% This function is designed to provide a utility for extracting object
% trajectory data from the raw DATA.globalTrajectories data set.

% Author: James A. Douthwaite

% Fetch states from Trajectory matrix
function [objectStateData,indexOfLastState] = OMAS_getTrajectoryData(trajectoryMatrix,globalIDvector,objectID,step)
% This is the generic function for returning the state trajectory data for
% a given agent ID number.
% INPUTS:
% trajectoryMatrix - a complete matrix of object global states [objects*n x steps].
% objectNumber     - The the index of desired objects ID in the SIM.OBJECTS list.
% step             - The requested simulation step.
%  + >steps        - Return states for all timesteps.
%  + nan           - Return states for all timesteps until the agent way idle.
%  + step          - Return a specific state at a given timestep.

% DETERMINE THE NUMBER OF OUTPUT STATES
numberOfGlobalStates = 10;                                                          % Determing the number of system states
systemStates = size(trajectoryMatrix,1);

% INPUT HANDLING
assert(mod(systemStates,numberOfGlobalStates) == 0,'The dimensions of the state matrix provided is incorrect.');
assert(any(globalIDvector == objectID),'The query objectID must belong to the global ID vector.');

% DEFINE THE INDEX OF THE OBJECT FROM ITS OBJECT ID 
logicalIndices = globalIDvector == objectID;                                        % The logical position of the ID
indexPosition = 0;
for i = 1:numel(logicalIndices)
   if logicalIndices(i)
      indexPosition = i; 
      break
   end
end

%INFER INDEX PROPERTIES
indexSet = uint16(zeros(2,systemStates/numberOfGlobalStates));
indexSet(1,:) = (0:numberOfGlobalStates:systemStates-numberOfGlobalStates) + 1;     % Declare the state indices
indexSet(2,:) = numberOfGlobalStates:numberOfGlobalStates:systemStates;

% DEFINE THE STATE LIMITS (INDICES) THAT 
objectStateIndices = indexSet(:,indexPosition);

% RETURN ALL TIMESERIES DATA
if step > size(trajectoryMatrix,2) 
    % FULL TIME SERIES
    objectStateData = trajectoryMatrix(objectStateIndices(1):objectStateIndices(2),:);   % Otherwise extract the complete timeseries
    indexOfLastState = uint32(size(objectStateData,2));                               % Indices of the last state
    return
end
% RETURN ALL STATES UPTO THE POINT OF INACTIVITY (INDICATED BY NaNs)
if isnan(step)
        % RETURN ONLY VALID STATES (ALL STATES UPTO TERMINATION 'NaN')
        objectStateData = trajectoryMatrix(objectStateIndices(1):objectStateIndices(2),:);  % Full timeseries
        % IF THE AGENT IS ACTIVE AT ALL TIMESTEPS
        logicalMatrix = logical(isnan(objectStateData(1:3,:)));
        if ~any(any(logicalMatrix))
            indexOfLastState = uint32(size(objectStateData,2));
            return
        end
        % THERE ARE NAN OCCURANCES
        indexOfLastState = size(logicalMatrix,2);
        for i = size(logicalMatrix,2):1
            if sum(logicalMatrix(:,i)) ~= 3
                indexOfLastState = size(logicalMatrix,2) - i;
                break
            end
        end
        objectStateData = objectStateData(:,1:indexOfLastState);                            % Final valid object state
        indexOfLastState = uint32(indexOfLastState);
else
    % SPECIFIC STEP
    objectStateData = trajectoryMatrix(objectStateIndices(1):objectStateIndices(2),step);% A specific time step has been specified
    indexOfLastState = uint32(step);
end

end
