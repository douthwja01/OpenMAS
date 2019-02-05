       % GET AGENT COMPUTATION STATISTICS 
       function [agentStatistics] = getAgentStatistics(obj,cycleData)
           % This function computes the agent statistics across all
           % Monte-Carlo cycles.
           
           % INPUT HANDLING
           if ~exist('cycleData','var')
              cycleData = obj.samples; 
           end
           if isempty(cycleData)
              warning('No agent data provided.');
              agentStatistics = [];
              return
           end
           
           fprintf('[%s]\tCalculating agent statistics...\n',obj.phase);
           
           % DATA CONTAINERS
           totalCycles = numel(cycleData);
           agentStatistics = struct('meanLoopTime',0,...
                                    'maxLoopTime',0,...
                                    'minLoopTime',0,...
                                    'meanLoopTimeSeries',[],...
                                    'meanTimeVector',[]);
           
           % CALCULATE THE EVENT STATISTICS
           for cycle = 1:totalCycles
               % // CALCULATE THE AGENT COMPUTATION MEANS FOR EACH CYCLE //
               metaObjects = cycleData(cycle).META.OBJECTS;
               agentMetaObjects = metaObjects([metaObjects.type] == OMAS_objectType.agent);              
               agentLogicals = ismember(cycleData(cycle).META.globalIDvector,[agentMetaObjects.objectID]);
               % GET THE AGENT-OBJECT SUBSET
               agentIndex  = cycleData(cycle).DATA.objectIndex(agentLogicals);
               agentNumber = cycleData(cycle).META.totalAgents;
               cycleTimeVector = cycleData(cycle).META.TIME.timeVector;
               
               % CYCLE DATA
               cycleMEANS = struct('algorithm_dt',[],...
                                   'algorithm_mean_dt',0,...
                                   'algorithm_max_dt',0,...
                                   'algorithm_min_dt',0);
                               
               % CALCULATE THE AGENT MEANS
               for agentNum = 1:agentNumber
                   % DATA LOCAL TO A GIVEN AGENT OF A GIVEN CYCLE
                   agentDATA = agentIndex{agentNum}.DATA;             
                   % THE MEAN NOMINAL COMPUTATION TIME
                   cycleMEANS.algorithm_mean_dt = cycleMEANS.algorithm_mean_dt + agentDATA.algorithm_mean_dt/agentNumber;
                   % THE MEAN MAXIMUM COMPUTATION TIME
                   cycleMEANS.algorithm_max_dt = cycleMEANS.algorithm_max_dt + agentDATA.algorithm_max_dt/agentNumber;
                   % THE MEAN MINIMUM COMPUTATION TIME
                   cycleMEANS.algorithm_min_dt = cycleMEANS.algorithm_min_dt + agentDATA.algorithm_min_dt/agentNumber;
                   
                   % BUILD CYCLE MEAN COMPUTATION TIME-TIMESERIES
                   % The simulations are of varying lengths, therefore the 
                   % computation times must be padded to the maximum 
                   % duration to ensure regular data.
                   algIndicator = zeros(size(cycleTimeVector));
                   algIndicator(1:numel(agentDATA.algorithm_indicator)) = agentDATA.algorithm_indicator;
                   algCompTime = zeros(size(cycleTimeVector));
                   algCompTime(1:numel(agentDATA.algorithm_indicator)) = agentDATA.algorithm_dt;
                   % USE THE PADDED TIME-SERIES TO AVERAGE THE MEANS.
                   computationTimeSeries = algCompTime;
                   
                   if isempty(cycleMEANS.algorithm_dt)
                       cycleMEANS.algorithm_dt = computationTimeSeries/agentNumber;
                   else
                       cycleMEANS.algorithm_dt = cycleMEANS.algorithm_dt + computationTimeSeries/agentNumber;
                   end
                   % LEGACY
%                    computationTimeSeries = agentDATA.algorithm_dt(agentDATA.algorithm_indicator == 1); % Timeseries data for only when the algorithm was ran
%                    computationTimeSeries = agentDATA.algorithm_dt;
%                    if isempty(cycleMEANS.algorithm_dt)
%                        cycleMEANS.algorithm_dt = computationTimeSeries/agentNumber;
%                    else
%                        cycleMEANS.algorithm_dt = cycleMEANS.algorithm_dt + computationTimeSeries/agentNumber;
%                    end
               end
               
               % /////////// CALCULATE THE CROSS-CYCLE MEANS //////////////
               % GLOBAL MEAN NOMINAL LOOP TIME
               agentStatistics.meanLoopTime = agentStatistics.meanLoopTime + cycleMEANS.algorithm_mean_dt/totalCycles;
               % GLOBAL MEAN MAXIMUM LOOP TIME
               agentStatistics.maxLoopTime = agentStatistics.maxLoopTime + cycleMEANS.algorithm_max_dt/totalCycles;
               % GLOBAL MEAN MINIMUM LOOP TIME
               agentStatistics.minLoopTime = agentStatistics.minLoopTime + cycleMEANS.algorithm_min_dt/totalCycles;
               
               % ///// BUILD CYCLE MEAN COMPUTATION TIME-TIMESERIES ///////
               computationTimeSeries = cycleMEANS.algorithm_dt;
               if isempty(agentStatistics.meanLoopTimeSeries)
                   % IF THIS IS THE FIRST CYCLE-SAMPLE
                   agentStatistics.meanLoopTimeSeries = computationTimeSeries/totalCycles;
               else
                   % IF THE CURRENT MEAN VECTOR IS LONGER THAN THE CYCLE-SAMPLE
                   % There is no gaurantee that the timeseries are of the
                   % sample length therefore we must pad the estimates.
%                    if numel(agentStatistics.meanLoopTimeSeries) > numel(computationTimeSeries)
%                        % CONTAINER OF THE LARGER SIZE
%                        newDataFrame = NaN(size(agentStatistics.meanLoopTimeSeries));
% %                        newDataFrame = zeros(size(agentStatistics.meanLoopTimeSeries));
%                        % EXTEND THE COMPUTATION TIMESERIES
%                        newDataFrame(1:numel(computationTimeSeries)) = computationTimeSeries;
%                        % MODIFY THE SIZE OF THE COMPUTATION TIME-SERIES
%                        computationTimeSeries = newDataFrame;
%                    else
%                        % CONTAINER OF THE LARGER SIZE
%                        newDataFrame = NaN(size(computationTimeSeries));
% %                        newDataFrame = zeros(size(computationTimeSeries));
%                        % EXTEND THE CURRENT SAMPLE MEAN
%                        newDataFrame(1:numel(agentStatistics.meanLoopTimeSeries)) = agentStatistics.meanLoopTimeSeries;
%                        % MODIFY THE SIZE OF THE MONTE-CARLO TIME-SERIES 
%                        agentStatistics.meanLoopTimeSeries = newDataFrame;
%                    end
                   % IF THE CURRENT MEAN VECTOR IS LONGER THAN THE CYCLE-SAMPLE
                   agentStatistics.meanLoopTimeSeries = agentStatistics.meanLoopTimeSeries + computationTimeSeries/totalCycles;
               end
           end
           
           % COMPUTE MEAN TIME VECTOR
           % The 'meanLoopTimeSeries' is now padded to ensure it is of the
           % same length as TIME.timeVector to yield better manipulation.
           agentStatistics.meanTimeVector = cycleData(end).DATA.timeVector;
           
%            agentStatistics.meanTimeVector = cycleData(end).DATA.timeVector(1:numel(agentStatistics.meanLoopTimeSeries));
           
           fprintf('[%s]\t... Complete.\n',obj.phase);
       end