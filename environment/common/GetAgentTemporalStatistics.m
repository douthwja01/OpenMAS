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