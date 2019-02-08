%% OPENMAS AGENT STATISTICS (OMAS_agentStatics.m) %%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates adds statistical data from the agent.DATA field
% to the output structure as part of the post-processing.

% Author: James A. Douthwaite 06/09/2017

function [objectIndex,MEANS] = OMAS_agentStatistics(SIM,objectIndex)
% Each event cant be identified by its .type and .name
% INPUTS:
% SIM         - Local copy of the simulation META variable
% objectIndex - 
% OUTPUTS:
% objectIndex - Updated 'objectIndex' with statistics

fprintf('[%s]\tParsing AGENT data...\n',SIM.phase);
% INPUT HANDLING
if isempty(objectIndex)
    fprintf('[%s]\tNo AGENT history to parse.\n',SIM.phase);
    return
end
% BUILD THE "MEANS" DATA STRUCTURE
MEANS = struct(); global_mean_dt = 0;
% Isolate the agent set
agentOBJECTS = objectIndex([SIM.OBJECTS(:).type] == OMAS_objectType.agent);  %(SIM.globalIDvector == objectIndex{ind}.objectID);
% MOVE THROUGH THE OBJECT SET 
for ind = 1:SIM.totalAgents
    % Check if a 'DATA' structure is present
    if isprop(agentOBJECTS{ind},'DATA') || isempty(objectIndex{ind}.DATA)
       continue 
    end
    % Query the agent.DATA property
    agentDATA = objectIndex{ind}.DATA; 
    
    % GET THE ALGORIHM TIMESERIES DATA & SUMMARY INFORMATION
    if isfield(agentDATA,'dt')
        [agentDATA] = get_algorithmTimeData(agentDATA);
        % STORE THE AGENT'S MEAN TIME
        fieldName = sprintf('objectID%d_mean_dt',objectIndex{ind}.objectID);
        MEANS.(fieldName) = agentDATA.mean_dt;
        % GET THE (GLOBAL) ALGORITHM COMPUTATION-TIME DATA
        global_mean_dt = global_mean_dt + agentDATA.mean_dt/SIM.totalAgents;
    end
    % Reassign the DATA to agent structure
    objectIndex{ind}.DATA = agentDATA;                                     % Update object-side data
end
end

% GENERATE THE TIMING PARAMETERS FROM THE ALGORITHM TIMESERIES
function [updatedAgentData] = get_algorithmTimeData(agentDATA)
% This function calculates all the timing parameters for the algorithm data
% deposited in the .DATA attribute of each agent.

% INITIALISE THE OUTPUT DATA CONTAINER
updatedAgentData = agentDATA;

if ~isfield(agentDATA,'indicator')
    % ASSUME ACTIVE CONSTANTLY
    agentDATA.indicator = ones(size(agentDATA.dt));
end

% DETERMINE THE ALGORITHM COMPUTATION TIME TIME-SERIES
executedSteps = logical(agentDATA.indicator);                              % Convert algorithm ran indicators into logicals
valid_algorithm_dt = agentDATA.dt(executedSteps);                          % Get the times where the computations were ran

% GET THE (AGENT) ALGORITHM COMPUTATION-TIME DATA
updatedAgentData.mean_dt = sum(valid_algorithm_dt)/numel(valid_algorithm_dt);
updatedAgentData.max_dt = max(valid_algorithm_dt);
updatedAgentData.min_dt = min(valid_algorithm_dt);
end