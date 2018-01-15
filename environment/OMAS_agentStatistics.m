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
% MOVE THROUGH THE OBJECT SET 
for ind = 1:numel(objectIndex)
    % GET THE OBJECTS META OBJECT
    METAObject = SIM.OBJECTS(SIM.globalIDvector == objectIndex{ind}.objectID); 
    isValidAgent = METAObject.type == OMAS_objectType.agent || 1 == isprop(objectIndex{ind},'DATA');
    if ~isValidAgent || isempty(objectIndex{ind}.DATA)
        continue
    end
    % QUERY THE AGENT.DATA PROPERTY
    agentDATA = objectIndex{ind}.DATA;
    % GET THE ALGORIHM TIMESERIES DATA & SUMMARY INFORMATION
    if isfield(agentDATA,'algorithm_dt')
        [agentDATA] = get_algorithmTimeData(agentDATA);
        % STORE THE AGENT'S MEAN TIME
        fieldName = sprintf('%s_mean',objectIndex{ind}.name);
        MEANS.(fieldName) = agentDATA.algorithm_mean_dt;
        % GET THE (GLOBAL) ALGORITHM COMPUTATION-TIME DATA
        global_mean_dt = global_mean_dt + agentDATA.algorithm_mean_dt/SIM.totalAgents;
    end
    objectIndex{ind}.DATA = agentDATA;                                     % Update object-side data
end

end

% GENERATE THE TIMING PARAMETERS FROM THE ALGORITHM TIMESERIES
function [updatedAgentData] = get_algorithmTimeData(agentDATA)
% This function calculates all the timing parameters for the algorithm data
% deposited in the .DATA attribute of each agent.

% INITIALISE THE OUTPUT DATA CONTAINER
updatedAgentData = agentDATA;

if ~isfield(agentDATA,'algorithm_indicator')
    % ASSUME ACTIVE CONSTANTLY
    agentDATA.algorithm_indicator = ones(size(agentDATA.algorithm_dt));
end

% DETERMINE THE ALGORITHM COMPUTATION TIME TIME-SERIES
executedSteps = logical(agentDATA.algorithm_indicator);                % Convert algorithm ran indicators into logicals
valid_algorithm_dt = agentDATA.algorithm_dt(executedSteps);            % Get the times where the computations were ran

% GET THE (AGENT) ALGORITHM COMPUTATION-TIME DATA
updatedAgentData.algorithm_mean_dt = sum(valid_algorithm_dt)/numel(valid_algorithm_dt);
updatedAgentData.algorithm_max_dt = max(valid_algorithm_dt);
updatedAgentData.algorithm_min_dt = min(valid_algorithm_dt);
end