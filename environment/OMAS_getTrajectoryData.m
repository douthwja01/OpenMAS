%% OPENMAS TRAJECTORY DATA EXTRACTION (simulation_getTrajectory.m) %%%%%%%%
% This function is designed to provide a utility for extracting object
% trajectory data from the raw DATA.globalTrajectories data set.

% Author: James A. Douthwaite

% Fetch states from Trajectory matrix
function [objectStateData,stateIndex] = OMAS_getTrajectoryData(DATA,objectID,step)
% This is the generic function for returning the state trajectory data for
% a given agent ID number.
% INPUTS:
% DATA     - The complete DATA structure
% objectID - The object ID to be indexed
% step     - The desired step
% OUTPUT:
% objectStateData - The object state evolution data

% DETERMINE THE NUMBER OF OUTPUT STATES
statecount = size(DATA.globalTrajectories,1)/DATA.totalObjects;            % Determing the number of system states

% IF THE INDEX IS NOT KNOWN
if ~isfield(DATA,'stateIndex')   
    systemStates = size(DATA.globalTrajectories,1);
    indexSet(1,:) = (0:statecount:systemStates-statecount) + 1;            % Declare the state indices
    indexSet(2,:) = statecount:statecount:systemStates;
    DATA.stateIndex = indexSet;
    % Clear redundant variables 
    clearvars indexSet systemStates
end
% DEFINE THE STATE LIMITS (INDICES) THAT 
stateLimits = DATA.stateIndex(:,objectID);

% EXTRACT DATA USING THE DATA.stateIndex
if exist('step','var') 
    if ischar(step) && strncmp(step,'last',4)
    	% LAST STATE REQUESTED
        objectStateData = DATA.globalTrajectories(stateLimits(1):stateLimits(2),:);   % Full timeseries
        validStates = ~any(isnan(objectStateData(1:3,:)));                            % non-NaN states
        stateIndex = find(validStates,1,'last');
        objectStateData = objectStateData(:,stateIndex);                              % Final valid object state
    elseif ischar(step) && strncmp(step,'valid',5)
        % ONLY VALID STATES (ALL STATES UPTO TERMINATION 'NaN')
        objectStateData = DATA.globalTrajectories(stateLimits(1):stateLimits(2),:);   % Full timeseries
        validStates = ~any(isnan(objectStateData(1:3,:)));                            % non-NaN states
        stateIndex = find(validStates,1,'last');
        objectStateData = objectStateData(:,1:stateIndex);                            % Final valid object state
    else
        % SPECIFIC STEP
        objectStateData = DATA.globalTrajectories(stateLimits(1):stateLimits(2),step);% A specific time step has been specified
        stateIndex = step;
    end
else
    % FULL TIME SERIES
    objectStateData = DATA.globalTrajectories(stateLimits(1):stateLimits(2),:);   % Otherwise extract the complete timeseries
    stateIndex = size(objectStateData,2);
end

clearvars newStateIndex stateIndex n statecount;
end

