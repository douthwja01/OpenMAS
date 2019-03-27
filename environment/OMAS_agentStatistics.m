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
MEANS = struct('global_mean_dt',0,...
               'global_max_dt',0,...
               'global_min_dt',0); 
           
% Isolate the agent set
agentOBJECTS = objectIndex([SIM.OBJECTS(:).type] == OMAS_objectType.agent);%(SIM.globalIDvector == objectIndex{ind}.objectID);
% MOVE THROUGH THE OBJECT SET 
for ind = 1:SIM.totalAgents
    
    % Check if a 'DATA' structure is present
    if ~isprop(agentOBJECTS{ind},'DATA') || isempty(agentOBJECTS{ind}.DATA)
       continue 
    end
    % Query the agent.DATA property
    agentDATA = agentOBJECTS{ind}.DATA; 
    
    % GET THE ALGORIHM TIMESERIES DATA & SUMMARY INFORMATION
    if isfield(agentDATA,'dt')
        % Get the analysis of the mean computation where the algorithm was
        % executed

        % Get the temporal parameters of the algorithm
        [mean_dt,max_dt,min_dt] = GetAgentTemporalStatistics(agentDATA.dt,agentDATA.indicator);
        
        % STORE THE AGENT'S MEAN TIME
        MEANS.(sprintf('ID%d_mean_dt',objectIndex{ind}.objectID)) = mean_dt;
        
        % GET THE (GLOBAL) ALGORITHM COMPUTATION-TIME DATA
        MEANS.global_mean_dt = MEANS.global_mean_dt + mean_dt/double(SIM.totalAgents);
        MEANS.global_max_dt  = MEANS.global_max_dt  + max_dt/double(SIM.totalAgents);
        MEANS.global_min_dt  = MEANS.global_min_dt  + min_dt/double(SIM.totalAgents);
%     else
%         warning('Agent %s[ID-%d] does not have a "dt" property for temporal analysis.',...
%             agentOBJECTS{ind}.name,agentOBJECTS{ind}.objectID);
    end
    % Reassign the DATA to agent structure
    objectIndex{[SIM.OBJECTS(:).objectID] == agentOBJECTS{ind}.objectID}.DATA = agentDATA;% Update object-side data                           
end
end

% GENERATE THE TIMING PARAMETERS FROM THE ALGORITHM TIMESERIES
function [mean_dt,max_dt,min_dt] = GetAgentTemporalStatistics(dt_vector,logical_vector)
% This function calculates all the timing parameters for the algorithm data
% deposited in the .DATA attribute of each agent.

% INPUTS:
% dt_vector     - The duration of computation at each time-step.
% action_vector - The logical vector indicating if the algorithm was ran.

% Input checking
assert(numel(dt_vector) == numel(logical_vector),'There must be an indicator logical for each time "dt" time-step.');

% DETERMINE THE ALGORITHM COMPUTATION TIME TIME-SERIES
valid_algorithm_dt = dt_vector(logical(logical_vector));                   % Get the times where the computations were ran

% GET THE (AGENT) ALGORITHM COMPUTATION-TIME DATA
mean_dt = sum(valid_algorithm_dt)/numel(valid_algorithm_dt);
max_dt  = max(valid_algorithm_dt);
min_dt  = min(valid_algorithm_dt);
end