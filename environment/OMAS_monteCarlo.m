%% OpenMAS MONTE-CARLO SIMULATION UTILITY (OMAS_monteCarlo.m) %%%%%%%%%%%%%
% This function is designed to compute a defined number of instances of the
% OpenMAS multi-agent simulator.
% TERMINOLOGY:
% study   - A set of initial conditions consiting of several sessions. 
% session - A defined set of initial conditions for which many OpenMAS 
%           instances are evaluated.
% cycle   - A individual instance of the OpenMAS scenario with the initial
%           conditions specified by the session settings.

classdef OMAS_monteCarlo
   properties
       % DEFAULT OMAS SIMULATION PARAMETERS
       OMAS_settings;                                                      % Container for the MC settings
       
       % DEFAULT MONTE-CARLO PARAMETERS 
       objects;                                % The object set
       maxObjectNumber = 300;                  % The maximum number of objects defining studies
       cycles;                                 % Default number of cycles
       sessions;                               % Parameter representing the num
       sessionLabels = {};                     % Parameter providing reference to each session
       phase;
       systemFiles = {'META','DATA','EVENTS','OBJECTS'};
       threadPool;                             % Complete in a asynchronous way
       shutDownOnComplete;                     % Shutdown on completion
       outputPath = pwd;
   end
   
   % ////////////////////////////// MAIN //////////////////////////////////
   methods
       % INITIALISER METHOD
       function [obj] = OMAS_monteCarlo(varargin)
           % INPUT HANDLING
           if length(varargin) == 1 && iscell(varargin)                    % Catch nested cell array inputs
               varargin = varargin{:};
           end
           
           % //////////////////// INPUT HANDLING ////////////////////////// 
           % PARSE MONTE-CARLO PARAMETERS
           [monteCarloConfig] = obj.configurationParser(obj.getDefaultMCconfig(),varargin);
           [obj] = obj.configurationParser(obj,monteCarloConfig);     
           % PARSE DEFAULT OMAS CONFIGURATION
           obj.OMAS_settings = obj.configurationParser(obj.getDefaultOMASConfig(),varargin); % Parse input config against the default config
           % CONFIRM PROPOSED OBJECTS ARE VALID
           assert(~isempty(obj.objects),'Please provide a object array of the form: [numberOfAlgorithms x numberOfObjects]');
           assert(numel(obj.objects) < obj.maxObjectNumber,sprintf('The number objects (%d) in the proposed study is to high.',numel(obj.objects))) 
           obj.sessions = size(obj.objects,1);
           % //////////////////////////////////////////////////////////////
           fprintf('\n[%s]\tIntialising Monte-Carlo instance.\n',obj.phase);
           % ENABLE MATLAB PAUSING
           pause on
       end
   end
   methods
       % EVALUATE THE PROPOSED CYCLE
       function [] = evaluateAllCycles(obj)
           % This function completes the proposed cycles given the settings
           % specified settings. This function assumes that we are 
           % evaluating:
           % VARIABLES:
           % - objectSets        - [m x n]
           % - numberPopulations - [a b c]
           % CONSTANTS:
           % - cycles
           % - SETTINGS
           
           % GENERATE THE MONTE-CARLO DIRECTORY
           [highestLevelPath] = obj.createMonteCarloDirectory(obj.outputPath);
           
           % CREATE SESSION LABELS
           for i = 1:obj.sessions
               sessionObjects = obj.objects(i,:);
               sessionObjects = sessionObjects(cellfun('isempty',sessionObjects) == 0); 
               
               % WRITE SESSION LABEL 
               % Simply take the class of the first object
               obj.sessionLabels{i} = class(sessionObjects{1});
           end
           
           % IF REQUIRED, ENSURE A POOL IS PRESENT
           if obj.threadPool
               obj.threadpoolInterface('CONFIRM');
           end
           
           % ////////////// POPULATE THE SESSION DIRECTORIES //////////////
           obj.phase = 'CYCLING';
           for i = 1:obj.sessions 
               % Each session represents a unique object set that will be
               % simulated over 'n' cycles.
               
               % //////////// DEFINE THE SESSION PARAMETERS ///////////////
               % OMIT EMPTY ENTRIES IN THE OBJECT INDEX
               sessionObjects = obj.objects(i,:);
               sessionObjects = sessionObjects(cellfun('isempty',sessionObjects) == 0);           
               sessionLabel = sprintf('-objects(%d)-label(%s)',numel(sessionObjects),obj.sessionLabels{i});
               % ALLOCATE NEW SESSION PATH FOR A GIVEN STUDY
               sessionPath = obj.createSessionDirectory(highestLevelPath,sessionLabel);                            
               % ///////////////// COMPUTE THE CYCLES /////////////////////
               if obj.threadPool
                   obj.parallelCycles(sessionPath,sessionObjects);
               else
                   obj.sequentialCycles(sessionPath,sessionObjects);
               end
               % //////////////////////////////////////////////////////////
               clearvars -except obj highestLevelPath
           end
           
           % //////////////////// DATA IMPORTATION ////////////////////////
           % The number of data points imported can potentially be very
           % large. For this reason, the parallel pool is retained.
           % The data that will be imported:
           obj.phase = 'ANALYSIS';
           
           % GET LIST OF OUTPUT DATA LOCATIONS BY SESSION
           [sessionList] = obj.getPathSubDirectories(highestLevelPath);
           sessionLabelVector = obj.sessionLabels;
           if obj.threadPool
               parfor i = 1:obj.sessions
                   % SESSION PATH
                   sessionPath = [highestLevelPath,'\',sessionList{i}];
                   % COMPUTE SESSION BASED STATISTICS 
                   summaryData(i) = obj.computeSessionAnalyitics(sessionPath);
                   % ADD GENERAL FIELDS
                   summaryData(i).label = sessionLabelVector{i};
               end
           else
               for i = 1:obj.sessions
                   % SESSION PATH
                   sessionPath = [highestLevelPath,'\',sessionList{i}];
                   % COMPUTE SESSION BASED STATISTICS 
                   summaryData(i) = obj.computeSessionAnalyitics(sessionPath);
                   % ADD GENERAL FIELDS
                   summaryData(i).label = sessionLabelVector{i};
               end
           end
           % /////// CALCULATE THE MONTE-CARLO INSTANCE STATISTICS ////////
           % Using the data collected from each of the monte-carlo sessions
           % we can calculate and compare statistics from all sessions.
           obj.phase = 'OUTPUT';
           
           studyOutputPath = [highestLevelPath,'\'];
           
           % GET THE OpenMAS FIGURE PROPERTIES STRUCTURE (for continuity)
           [figureProperties] = OMAS_figureProperties();
                     
           % GENERATE THE MEAN COMPUTATIONAL TIMES
           obj.get_meanComputationTimeSeries(figureProperties,studyOutputPath,summaryData);
           % GENERATE THE FIGURES COMPARING MEAN COMPUTATIONAL TIMES TO AGENT POPULATION
           obj.get_meanComputationTimesVsPopulation(figureProperties,studyOutputPath,summaryData);
           % GENERATE THE FIGURES COMPARING MEAN NUMBER OF COLLISIONS TO AGENT POPULATION
           obj.get_meanCollisionsVsPopulation(figureProperties,studyOutputPath,summaryData);
           % //////////////////////////////////////////////////////////////
           
           % KILL THE POOL WHEN ALL SESSIONS/STUDIES ARE COMPLETE
           if obj.threadPool
               obj.threadpoolInterface('KILL');
           end
       end
       % COMPUTE ANALYSIS OVER A SESSIONS - CYCLE SET
       function [summaryData] = computeSessionAnalyitics(obj,sessionPath)
           % This function executes a series of basic analyses based on the
           % populated cycle directories. This is essentially the mean data
           % from the cycles described by the session conditions.
           
           % GET THE VALID SESSIONS CREATED
           [cycleDirectories] = obj.getValidCycleDirectories(sessionPath);
           totalSamples = numel(cycleDirectories);
           
           % CREATE THE SUMMARY CSV
           fileName = [sessionPath,'\','summaryData.csv'];
           headerList = {'collisions','waypoints','mean_dt','min_dt','max_dt','minSeparation','maxSeparation'};
           obj.pushDataToCSV(fileName,1,headerList)
           
           % CREATE CONTAINER
           summaryData = struct(...
               'label',[],...
               'cycles',totalSamples,...
               'agents',[],...
               'waypoints',[],...
               'mean_collisions',0,...
               'mean_waypoints',0,...
               'mean_minSeparation',0,...
               'mean_maxSeparation',0,...
               'mean_loopTime',0,...
               'mean_minLoopTime',0,...
               'mean_maxLoopTime',0,...
               'mean_timeVector',0,...
               'dt_timeSeries',0);
           
           % GET THE DATA FROM THE FIRST CYCLE
           [cycleSample,~] = obj.importPathData(obj.systemFiles,[sessionPath,'\',cycleDirectories{1}]);          
           % GET THE NUMBER OF AGENTS IN THE SIMULATION
           % This is only valid if all cycles are computed with the same
           % OMAS_settings, which they should be.
           summaryData.agents = sum([cycleSample.META.OBJECTS(:).type] == OMAS_objectType.agent);
           summaryData.waypoints = sum([cycleSample.META.OBJECTS(:).type] == OMAS_objectType.waypoint);
           
           clearvars cycleSample
            
           % /////////////// BEGIN RUNNING THROUGH CYCLES /////////////////
           for i = 1:totalSamples
               % THE PATH CYCLES
               cyclePath = [sessionPath,'\',cycleDirectories{i}];
               % IMPORT THE CYCLE SAMPLE
               [cycleSample,~] = obj.importPathData(obj.systemFiles,cyclePath);
               
               % //////////////// COMPUTE CYCLE DATA //////////////////////
               % COMPUTE THE AGENT-TIMING STATISTICS
               [mean_dt,min_dt,max_dt,dt_timeSeries,timeVector] = obj.cycle_timingAnalytics(cycleSample.META,cycleSample.OBJECTS);
               % COMPUTE THE AGENT-TRAJECTORY STATISTICS
               [minSeparation,maxSeparation] = obj.cycle_trajectoryAnalytics(cycleSample.META,cycleSample.DATA);
               % COMPUTE THE CYCLE - EVENT BASED STATISTICS 
               [collisions,waypoints] = obj.cycle_eventAnalytics(cycleSample.DATA);
               % //////////////////////////////////////////////////////////
               
               % PUSH CYCLE DATA TO .CSV (for reference)
               obj.pushDataToCSV(fileName,0,[collisions,waypoints,mean_dt,min_dt,max_dt,minSeparation,maxSeparation])
               
               % //////////////// UPDATE MEAN VALUES //////////////////////
               summaryData.mean_collisions     = summaryData.mean_collisions + collisions/totalSamples;
               summaryData.mean_waypoints      = summaryData.mean_waypoints + waypoints/totalSamples;
               summaryData.mean_minSeparation  = summaryData.mean_minSeparation + minSeparation/totalSamples;
               summaryData.mean_maxSeparation  =  summaryData.mean_maxSeparation + maxSeparation/totalSamples;
               % COMPUTE THE MEAN TIMING PARAMETERS
               summaryData.mean_loopTime    = summaryData.mean_loopTime + mean_dt/totalSamples;
               summaryData.mean_minLoopTime = summaryData.mean_minLoopTime + min_dt/totalSamples;
               summaryData.mean_maxLoopTime = summaryData.mean_maxLoopTime + max_dt/totalSamples;
               summaryData.mean_timeVector  = timeVector;
               summaryData.dt_timeSeries    = summaryData.dt_timeSeries + dt_timeSeries/totalSamples;
               
               % CLEAR REDUNDANT VARIABLES
               clearvars -except obj sessionPath cycleDirectories totalSamples fileName summaryData
           end
           % //////////// PUSH SESSION MEAN DATA TO .mat FILE /////////////
           save([sessionPath,'\meanData.mat'], '-struct', 'summaryData');
           clearvars -except summaryData;
       end
   end
   
    %% ////////////////////////// ANALYTICS /////////////////////////////// 
   methods (Static)
       % COMPUTE THE CYCLE AGENT TRAJECTORY STATISTICS
       function [minSeparation,maxSeparation] = cycle_trajectoryAnalytics(cycle_META,cycle_DATA)
           % This function is designed to recover relevant trajectory
           % statistics and return them as representative samples.

           % BUILD OBJECT SET FOR ALL NON-WAYPOINTS
           collidableMETA = cycle_META.OBJECTS([cycle_META.OBJECTS.type] ~= OMAS_objectType.waypoint);
           
           [LIA,~] = ismember(cycle_META.globalIDvector,[collidableMETA.objectID]);
           collidableIDs = cycle_META.globalIDvector(LIA);
           
           % OUTPUT CONTAINER
           separationTimeSeries = inf(numel(collidableIDs),numel(collidableIDs),cycle_META.TIME.numSteps);
           for IDnumA = 1:numel(collidableIDs)
               % THE AGENT STATES
               objectStatesA = OMAS_getTrajectoryData(...
                   cycle_DATA.globalTrajectories,...
                   cycle_META.globalIDvector,...
                   collidableIDs(IDnumA),...
                   inf);

               for IDnumB = 1:numel(collidableIDs)
                   % OMIT SEPERATIONS BETWEEN ITSELF
                   if IDnumA ~= IDnumB 
                       % GET THE AGENT STATE TIMESERIES
                       objectStatesB = OMAS_getTrajectoryData(...
                           cycle_DATA.globalTrajectories,...
                           cycle_META.globalIDvector,...
                           collidableIDs(IDnumB),...
                           inf);
                       centroidSeparations = objectStatesB(1:3,:) - objectStatesA(1:3,:);  % seperation of the centroids
                       centroidSeparations = sqrt(sum(centroidSeparations.^2,1));
                       % STORE IN SEPERATION TIMESERIES MATRIX
                       separationTimeSeries(IDnumB,IDnumA,:) = centroidSeparations;
                   end
               end
           end
           % MIN AND MAXIMUM SEPERATION MATRICES [IDA by IDB]
           minABMatrix = min(separationTimeSeries,[],3);                   % Minimum agent-object separations over the timeseries
           maxABMatrix = max(separationTimeSeries,[],3);                   % Maximum agent-object separations over the timeseries
           % GET THE MINIMUM SEPARATION FOR ALL COLLIDABLES FOR THE COMPLETE TIMESERIES
           minSeparation = min(minABMatrix(~isinf(minABMatrix)));
           maxSeparation = max(maxABMatrix(~isinf(maxABMatrix)));
       end
       % COMPUTE THE CYCLE AGENT TIMING STATISTICS
       function [mean_dt,min_dt,max_dt,dt_timeSeries,timeVector] = cycle_timingAnalytics(cycle_META,cycle_OBJECTS)
           % GET AGENT COMPUTATION STATISTICS
           % This function computes the agent statistics across all
           % Monte-Carlo cycles.
           
           % // CALCULATE THE AGENT COMPUTATION MEANS FOR EACH CYCLE //
%            metaObjects = cycle_META.OBJECTS;
%            agentMetaObjects = metaObjects([metaObjects.type] == OMAS_objectType.agent);          
          
           % GET THE AGENT-OBJECT SUBSET
           agentNumber = double(cycle_META.totalAgents);
           timeVector  = cycle_META.TIME.timeVector;
           agentIndex  = cycle_OBJECTS([cycle_META.OBJECTS(:).type] == OMAS_objectType.agent);      
           
           % OUTPUTS
           mean_dt = 0;
           max_dt  = 0;
           min_dt  = 0;
           dt_timeSeries = zeros(size(timeVector));
           
           % CALCULATE THE AGENT MEANS
           for agentNum = 1:agentNumber
               % DATA LOCAL TO A GIVEN AGENT OF A GIVEN CYCLE
               agentDATA = agentIndex{agentNum}.DATA;
               % Agent temporal statistics
               [agent_mean_dt,agent_max_dt,agent_min_dt] = GetAgentTemporalStatistics(...
                   agentDATA.dt,...
                   agentDATA.indicator);             
               
               % THE MEAN NOMINAL COMPUTATION TIME
               mean_dt = mean_dt + agent_mean_dt/agentNumber;
               % THE MEAN MAXIMUM COMPUTATION TIME
               max_dt  = max_dt + agent_max_dt/agentNumber;
               % THE MEAN MINIMUM COMPUTATION TIME
               min_dt  = min_dt + agent_min_dt/agentNumber;
               
               % BUILD CYCLE MEAN COMPUTATION TIME-TIMESERIES
               % The simulations are of varying lengths, therefore the
               % computation times must be padded to the maximum
               % duration to ensure regular data.
               % PAD THE ALGORITHM TIMERS TO THE FULL TIMEVECTOR LENGTH
               algCompTime = zeros(size(timeVector));
               algCompTime(1:numel(agentDATA.indicator)) = agentDATA.dt;
               
               % USE THE PADDED TIME-SERIES TO AVERAGE THE MEANS.
               computationTimeSeries = algCompTime;
               if isempty(dt_timeSeries)
                   dt_timeSeries = computationTimeSeries/agentNumber;
               else
                   dt_timeSeries = dt_timeSeries + computationTimeSeries/agentNumber;
               end
           end
           
           % COMPUTE MEAN TIME VECTOR
           % The 'meanLoopTimeSeries' is now padded to ensure it is of the
           % same length as TIME.timeVector to yield better manipulation.
       end
       % COMPUTE THE CYCLE EVENT STATISTICS
       function [collisions,waypoints] = cycle_eventAnalytics(cycle_DATA)
           % This function extracts the 'EVENT' based data from a completed
           % OpenMAS cycle. 
           
           assert(isstruct(cycle_DATA),'The sample does not have a DATA sub-structure.');
                          
           % GET THE NUMBER OF CYCLE COLLISIONS
           collisions = double(cycle_DATA.collisions);
           % GET THE NUMBER OF CYCLE WAY-POINTS ACHIEVED
           waypoints  = double(cycle_DATA.waypointsAchieved);
       end
   end
   %% ///////////////////// DATA HANDLING UTILITIES ///////////////////////
   methods (Static)
       % ////////////////////// FIGURE GENERATION /////////////////////////
       % GET THE COMPARISON OF COLLISIONS TO POPULATION
       function [figureHandle] = get_meanCollisionsVsPopulation(figureProperties,figurePath,meanData)
           % This function computes the mean number of collisions against
           % on population.
                      
           % CONFIGURE THE PLOT ATTRIBUTES
           figurePath = strcat(figurePath,'mean_collisionsVspopulation');
           figureHandle = figure('Name','OpenMAS Collisions vs Population');
           set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
           set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
           ax = axes(figureHandle);
           
           % REORDER THE DATA BY NUMBER OF AGENTS
           [~,orderIndices] = sort([meanData.agents]);
           meanData = meanData(orderIndices);
           invalidSessionLabels = {};
           for i = 1:numel(meanData)
               % IF THE LABEL HAS ALREADY BEEN COMPARED
               if any(ismember(invalidSessionLabels,meanData(i).label))
                   continue
               end
               % ADD THE LABEL TO PREVIOUSLY COMPARED LIST
               invalidSessionLabels = horzcat(invalidSessionLabels,meanData(i).label);
               % /////////////// COMPARATIVE MEAN DATAS ///////////////////
               % THE LOGICALS FOR THAT LABEL
               labelLogical = ismember({meanData.label},meanData(i).label);
               labelSeries = meanData(labelLogical);
               % PREPARE THE SERIES COMPUTATION TIMES SERIES             
               collisionSeries = [labelSeries(:).mean_collisions];
               populationSeries = [labelSeries(:).agents];
               % Define label
               labelString = ['$',strrep(labelSeries(1).label,'_','-'),'$'];
               % DEFINE ERROR BAR PLOT
               hold on;
               plot(ax,populationSeries,collisionSeries,...
                        'LineWidth',figureProperties.lineWidth,...  
                        'Color',rand(3,1),...
                        'DisplayName',labelString,...
                        'Marker','o',...
                        'MarkerSize',10);
           end
           % FIGURE PROPERTIES
           grid on; box on; grid minor;
           % Title
           title(ax,...
               'Mean Collisions vs Population',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontsize',figureProperties.titleFontSize,...
               'FontSmoothing','on');
           % X-Label
           xlabel(ax,...
               'Agent population (n)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           % Y-Label
           ylabel(ax,...
               'Mean collisions (n)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           % Legend
           le = legend(ax,'Location','northeast');
           set(le,'Interpreter',figureProperties.interpreter);
           % Axes
           set(ax,...
               'TickLabelInterpreter',figureProperties.interpreter,...
               'Fontweight',figureProperties.fontWeight,...
               'FontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on',...
               'Color',figureProperties.axesColor,...
               'GridLineStyle','--',...
               'GridAlpha',0.25,...
               'GridColor','k');
           % Cycle annotation
           cycleHandle = annotation('textbox',[0.025 0.025 0.15 0.04],...
               'String',sprintf('Cycles: %d',meanData(1).cycles),...
               'FitBoxToText','off');
           set(cycleHandle,...
               'interpreter',figureProperties.interpreter,...
               'FontSize',figureProperties.axisFontSize);
           
           hold off;
           % SAVE THE OUTPUT FIGURE
           savefig(figureHandle,figurePath);
           % PUBLISH TO PDF
           if figureProperties.publish
               set(figureHandle,'Units','Inches');
               pos = get(figureHandle,'Position');
               set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
               print(figureHandle,figurePath,'-dpdf','-r0');
           end
       end
       % GET THE COMPARISON Of MEAN COMPUTATION TIME AND AGENT NUMBER
       function [figureHandle] = get_meanComputationTimesVsPopulation(figureProperties,figurePath,meanData)
           % This function computes the mean computation times of each
           % session against the number of agents in its population.
                      
           % CONFIGURE THE PLOT ATTRIBUTES
           figurePath = strcat(figurePath,'mean_computationTimeVspopulation');
           figureHandle = figure('Name','OpenMAS Computation Time vs Population');
           set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
           set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
           ax = axes(figureHandle);
           
           % REORDER THE DATA BY NUMBER OF AGENTS
           [~,orderIndices] = sort([meanData.agents]);
           meanData = meanData(orderIndices);
           invalidSessionLabels = {};
           for i = 1:numel(meanData)
               % IF THE LABEL HAS ALREADY BEEN COMPARED
               if any(ismember(invalidSessionLabels,meanData(i).label))
                   continue
               end
               % ADD THE LABEL TO PREVIOUSLY COMPARED LIST
               invalidSessionLabels = horzcat(invalidSessionLabels,meanData(i).label);
               % /////////////// COMPARATIVE MEAN DATAS ///////////////////
               % THE LOGICALS FOR THAT LABEL
               labelLogical = ismember({meanData.label},meanData(i).label);
               labelSeries = meanData(labelLogical);
               % PREPARE THE SERIES COMPUTATION TIMES SERIES             
               minSeries  = [labelSeries(:).mean_minLoopTime];
               maxSeries  = [labelSeries(:).mean_maxLoopTime];
               meanSeries = [labelSeries(:).mean_loopTime];
               populationSeries = [labelSeries(:).agents];
               % DEFINE ERROR BAR PLOT
               hold on;
               labelString = ['$',strrep(labelSeries(1).label,'_','-'),'$'];
               errorbar(ax,populationSeries,meanSeries,minSeries,maxSeries,...
                        'LineWidth',figureProperties.lineWidth,...  
                        'Color',rand(3,1),...
                        'DisplayName',labelString,...
                        'Marker','o',...
                        'MarkerSize',10);
           end
           % FIGURE PROPERTIES
           grid on; box on; grid minor;
           % Title
           title(ax,...
               'Mean Computational Times vs Population',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontsize',figureProperties.titleFontSize,...
               'FontSmoothing','on');
           % X-Label
           xlabel(ax,...
               'Agent population (n)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           % Y-Label
           ylabel(ax,...
               'Computation Time (s)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           % Legend
           le = legend('Location','northeast');
           set(le,'Interpreter',figureProperties.interpreter);
           % Axes
           set(ax,...
               'TickLabelInterpreter',figureProperties.interpreter,...
               'Fontweight',figureProperties.fontWeight,...
               'FontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on',...
               'Color',figureProperties.axesColor,...
               'GridLineStyle','--',...
               'GridAlpha',0.25,...
               'GridColor','k');
           % Cycle annotation
           cycleHandle = annotation('textbox',[0.025 0.025 0.15 0.04],...
               'String',sprintf('Cycles: %d',meanData(1).cycles),...
               'FitBoxToText','off');
           set(cycleHandle,...
               'interpreter',figureProperties.interpreter,...
               'FontSize',figureProperties.axisFontSize);
           
           hold off;
           % SAVE THE OUTPUT FIGURE
           savefig(figureHandle,figurePath);
           % PUBLISH TO PDF
           if figureProperties.publish
               set(figureHandle,'Units','Inches');
               pos = get(figureHandle,'Position');
               set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
               print(figureHandle,figurePath,'-dpdf','-r0');
           end
       end
       % GET THE COMPARISON OF MEAN COMPUTATION TIME-SERIES DATA
       function [figureHandle] = get_meanComputationTimeSeries(figureProperties,figurePath,meanData)
           % This function computes the mean computation time-series of all
           % agents and compares them with each session.
           
           % ASSUMPTIONS:
           % meanData - The mean times are given in seconds.
           
           
           % CONFIGURE THE PLOT ATTRIBUTES
           figurePath = strcat(figurePath,'mean_computationTimes');
           figureHandle = figure('Name','OpenMAS Mean Computational Timeseries');
           set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
           set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
           ax = axes(figureHandle);
                      
           totalSamples = numel(meanData);
           for sessionNum = 1:totalSamples
               hold on;
               % NAME BASED ON SESSION
               displayName = ['n=',num2str(meanData(sessionNum).agents),',$ $',meanData(sessionNum).label];
               displayName = ['$',strrep(displayName,'_','-'),'$'];
               % PLOT THE MEAN TRACES FOR EACH SESSION
               plot(ax,...
                    meanData(sessionNum).mean_timeVector',...
                    meanData(sessionNum).dt_timeSeries',...
                    'Color',rand(3,1),...
                    'LineWidth',figureProperties.lineWidth,...
                    'DisplayName',displayName);
           end          
           grid on; box on; grid minor;
           title(ax,...
               'Mean Computational Timeseries',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.titleFontSize,...
               'FontSmoothing','on');
           xlabel(ax,...
               't (s)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           ylabel(ax,...
               'Computation Time (s)',...
               'interpreter',figureProperties.interpreter,...
               'fontweight',figureProperties.fontWeight,...
               'fontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on');
           % Legend
           le = legend(ax,'Location','northeast');
           set(le,'Interpreter',figureProperties.interpreter);
           % Axes
           set(ax,...
               'TickLabelInterpreter',figureProperties.interpreter,...
               'Fontweight',figureProperties.fontWeight,...
               'FontSize',figureProperties.axisFontSize,...
               'FontSmoothing','on',...
               'Color',figureProperties.axesColor,...
               'GridLineStyle','--',...
               'GridAlpha',0.25,...
               'GridColor','k');
           % Cycle annotation
           cycleHandle = annotation('textbox',[0.025 0.025 0.15 0.04],...
               'String',sprintf('Cycles: %d',meanData(1).cycles),...
               'FitBoxToText','off');
           set(cycleHandle,...
               'interpreter',figureProperties.interpreter,...
               'FontSize',figureProperties.axisFontSize);        
           hold off;
           % SAVE THE OUTPUT FIGURE
           savefig(figureHandle,figurePath);
           % PUBLISH TO PDF
           if figureProperties.publish
               set(figureHandle,'Units','Inches');
               pos = get(figureHandle,'Position');
               set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
               print(figureHandle,figurePath,'-dpdf','-r0');
           end
       end
       
       % ////////// EXPORTING //////////
       % PUSH VALUES TO CSV FILE
       function pushDataToCSV(filename,isHeaderLogical,rowData)
           % This function is designed to provide a mechanism for entering
           % values into a csv file with a designated header on the first
           % row.
           
           if isHeaderLogical && iscell(rowData)
                % //////// write the header string to the file ////////////
                %turn the headers into a single comma seperated string if it is a cell
                %array, 
                r = 0;
                header_string = rowData{1};
                for i = 2:length(rowData)
                    header_string = [header_string,',',rowData{i}];
                end
                %if the data has an offset shifting it right then blank commas must
                %be inserted to match
                if r>0
                    for i=1:r
                        header_string = [',',header_string];
                    end
                end

                %write the string to a file
                fid = fopen(filename,'w');
                fprintf(fid,'%s\r\n',header_string);
                fclose(fid);
           else
               % APPEND THE HEADER TO THE FILE
%                fid = fopen(filename,'a');
               % PUSH MESSAGE TO FILE
%                fprintf(fid,rowData);
               c = 0;
               r = 0;
               dlmwrite(filename, rowData,'-append','delimiter',',','roffset', r,'coffset',c);
           
           end
%                 cycleData = [collisions,waypoints,mean_dt,min_dt,max_dt];
%                dlmwrite(fileName,cycleData,'delimiter',',','-append');
           
       end
       % ////////// IMPORTING //////////
       % IMPORT OMAS FILES FROM A SPECIFIC SESSION DIRECTORY
       function [cycleSample,isSuccessful] = importPathData(OMASFileNames,cyclePath)
           % This function imports the specified ".mat" files from the
           % specified path. This is typically the output path of a given
           % OpenMAS cycle.
           
           cycleSample = [];
           isSuccessful = 0;
           
           % INPUT HANDLING
           if iscell(cyclePath)
               cyclePath = char(cyclePath);
           elseif ~ischar(cyclePath)
               warning('Please provide a valid path to the session directory');
               return
           end

           % IMPORT THE SESSION .MAT FILES
           cycleSample = cell2struct(cell(size(OMASFileNames)),OMASFileNames,2);
           
           % FOR EACH KNOWN DATA OBJECT
           for fileName = 1:numel(OMASFileNames)
               % PATH TO OMAS VARIABLES
               filePath = strcat(cyclePath,'\',OMASFileNames{fileName},'.mat');
               try
                   load(filePath,'-mat');
                   % IMPORT FILE AS VARIABLE
                   cycleSample.(OMASFileNames{fileName}) = eval(OMASFileNames{fileName});
               catch importError
                   warning('Unable to import data from:%s',filePath);
                   warning(importError.message);
               end
           end
           
           % DATA HAS BEEN SUCCESSFULLY IMPORTED
           if ~isempty(cycleSample)
               isSuccessful = 1;
           end
           
           % SEPERATE THE DATA STREAMS
           clearvars -except cycleSample isSuccessful
       end
       % GET LIST OF ALL VALID SESSION DATAS IN A DIRECTORY
       function [directorySet]    = getValidCycleDirectories(pathString)
           % This function is used to infer an analysis parameters given
           % the files in a provided path. It is assumed that the provided
           % path contains an index of session data sets.
           
           % INPUT CHECKING
           assert(ischar(pathString),'Please provide a valid path string.');
           
           % FETCH SESSION FOLDER SET
           rawDirectories = dir(pathString);
           % OMIT KNOWN BAD PATHS
           directorySet = rawDirectories(~ismember({rawDirectories.name},{'.','..'}));
           % OMIT NON-FOLDERS ( JUST THE LIST OF SESSION DIRECTORIES )
           directorySet = directorySet([directorySet.isdir] == 1); 
           % GET THE PATH STRINGS
           % GET ONLY CORRECTLY LABELLED SESSION DIRECTORIES
           directorySet = directorySet(strncmpi({directorySet.name},'cycledata',7));
           % CONVERT TO lIST OF PATHS
           directorySet = {directorySet.name}';
       end
       % GET LIST OF SUB-DIRECTORIES
       function [pathDirectories] = getPathSubDirectories(pathString)
           % INPUT CHECK
           assert(ischar(pathString),'Please provide a valid path string.');
           % FETCH SESSION FOLDER SET
           rawDirectories = dir(pathString);
           % OMIT KNOWN BAD PATHS
           pathDirectories = rawDirectories(~ismember({rawDirectories.name},{'.','..'}));
           % OMIT NON-FOLDERS ( JUST THE LIST OF SESSION DIRECTORIES )
           pathDirectories = pathDirectories([pathDirectories.isdir] == 1); 
           % GET THE PATH STRINGS
           pathDirectories = {pathDirectories.name}';
           % CLEAR VARIABLES
           clearvars -except pathDirectories
       end
   end
   
   %% ///////////////////////// CYCLE FUNCTIONS ///////////////////////////
   methods
       % RUN A DEFINED NUMBER OF SEQUENCIAL CYCLES
       function [obj] = sequentialCycles(obj,sessionPath,objectIndex)
           % This function computes a monte-carlo sequence of OMAS
           % instances.
                     
           % OPEN LOG (Only one handle is needed)
           logPath = strcat(sessionPath,'\log.txt'); 
           
           % /////////////// BEGIN THE CYCLE ITERATIONS ///////////////////
           fprintf('\n[%s]\tEXECUTING %d SEQUENTIAL MONTE-CARLO CYCLES...\n\n',obj.phase,obj.cycles);
           for cycle = 1:obj.cycles
               fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',obj.phase,cycle);
               [cycleMessage] = obj.singleCycle(sessionPath,obj.OMAS_settings,objectIndex);
               fprintf('\n[%s]\tCYCLE %d COMPLETE\n\n',obj.phase,cycle);          
               % MESSAGE
               message = sprintf('(Cycle:%d)\t%s',cycle,cycleMessage);
               % PUSH CYCLE SUMMARY TO LOG FILE
               obj.pushMessageToLog(logPath,message)
           end
           fprintf('[%s]\t%d SEQUENTIAL CYCLES COMPLETE.\n',obj.phase,obj.cycles);
           % //////////////////////////////////////////////////////////////
                      
           % ENSURE CLEARED
           clearvars -except obj
       end
       % RUN A DEFINED NUMBER OF PARALLEL CYCLES
       function [obj] = parallelCycles(obj,sessionPath,objectIndex)
           % Compute the proposed number of cycles in a parallel manner.
           % The resulting data generated by each simulation instance is
           % output to each session directory.
           
           fprintf('[%s]\tConfirming parallel worker pool...\n',obj.phase);
           
           % CHECK THE WORKER POOL IS INITIALISED
           obj.threadpoolInterface('CONFIRM');                             % Confirm the pool
            
           % OPEN LOG (Only one handle is needed)
           logPath = strcat(sessionPath,'\log.txt'); 
           
           % /////////////// BEGIN THE CYCLE ITERATIONS ///////////////////
           fprintf('\n[%s]\tEXECUTING %s PARALLEL MONTE-CARLO CYCLES...\n\n',obj.phase,num2str(obj.cycles));
           parfor cycle = 1:obj.cycles
               fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',obj.phase,cycle);
               [cycleMessage] = obj.singleCycle(sessionPath,obj.OMAS_settings,objectIndex);
               fprintf('\n[%s]\tCYCLE %d COMPLETE.\n\n',obj.phase,cycle);
               % MESSAGE
               message = sprintf('(Cycle:%d)\t%s',cycle,cycleMessage);
               % PUSH CYCLE SUMMARY TO LOG FILE
               obj.pushMessageToLog(logPath,message)
           end
           fprintf('[%s]\t%d PARALLEL CYCLES COMPLETE.\n',obj.phase,obj.cycles);
           % //////////////////////////////////////////////////////////////
                                 
           % ENSURE CLEARED
           clearvars -except obj
       end
   end
   methods (Static)
       % RUN ONE PARAMETERISED SIMULATION CYCLE
       function [cycleMessage] = singleCycle(sessionPath,SETTINGS,objectIndex)
           % This function computes a single monte-carlo simulation cycle
           % using the monte-carlo defined simulation parameters to create 
           % an instance of the OpenMAS simulation.
           % INPUTS:
           % sessionPath - The output directory of the cycles data.
           % SETTINGS    - The OpenMAS configurations.
           % objectIndex - The objects defining the scenario.
           
           % //////////////// PROCESS THE OBJECT INDEX ////////////////////  
           % PERTURB THE INITIAL XY POSITIONS (VALID FOR 2D & 3D)
           for objectNum = 1:numel(objectIndex)
               try               
                   % PERTURB THE GLOBAL POSITIONS
                   objectIndex{objectNum}.VIRTUAL.globalPosition = ...
                   + objectIndex{objectNum}.VIRTUAL.globalPosition ...
                   + [SETTINGS.positionSigma*randn(2,1);0];                % 2D perturbation;
               catch
                   warning('Unable to perturb global position of object %s',num2str(objectNum));
               end
           end
           clearvars -except sessionPath SETTINGS objectIndex
           % //////////////////////////////////////////////////////////////
           
           % ////// CREATE OpenMAS INSTANCE WITH THE DEFINED SETTINGS /////
           try
               OMAS_initialise(...
                   'duration',SETTINGS.duration,...
                   'dt',SETTINGS.dt,...                                        % Timing parameters
                   'figures',SETTINGS.figures,...
                   'warningDistance',SETTINGS.warningDistance,...
                   'conditionTolerance',SETTINGS.conditionTolerance,...
                   'visabilityModifier',SETTINGS.visabilityModifier,...
                   'threadPool',SETTINGS.threadPool,...                        % Disable the agent-loop thread pool
                   'verbosity',SETTINGS.verbosity,...
                   'gui',SETTINGS.gui,...
                   'monteCarloMode',SETTINGS.monteCarloMode,...
                   'outputPath',sessionPath,...
                   'objects',objectIndex);                                     % Enforce monte-carlo file creation
               % COMPLETION MESSAGE
               cycleMessage = 'Cycle completed successfully.';
           catch cycleError
               % FAILURE MESSAGE
               warning('\n...cycle error.\n');
               cycleMessage = cycleError.message;
           end
           clearvars -except cycleMessage;
           % //////////////////////////////////////////////////////////////
           
           % PAUSE FOR STABILITY
           pause(0.1); % 0.1s delay
       end    
   end

   %% //////////////////// MONTE-CARLO SETUP UTILITIES ////////////////////
   methods (Static)
       % ////////////////////// FILE OPERATIONS ///////////////////////////
       % PUSH MESSAGE TO LOG
       function pushMessageToLog(sessionPath,message)
           % This function takes a string resulting from a session and
           % generates the corresponding entry in the specified file.
           
           % OPEN LOG (Only one handle is needed)
           log = fopen(sessionPath, 'a');
           
           timeStamp = datestr(datetime('now'),'HH-MM-SS');
           entryString = sprintf('[%s]\t%s\n',timeStamp,char(message));
           % PUSH MESSAGE TO FILE
           fprintf(log,entryString);

           % CLOSE LOG
           fclose(log);     % For compatability 
       end
       % BUILD THE RELATIVE SESSION DIRECTORY
       function [sessionPath]    = createSessionDirectory(outputPath,fileLabel)
           % This function build the monte-carlo data output directory. The
           % naming of the file is of the form:
           % 'sessiondata 2018-09-28 Aobjects Bcycles
           
           % INPUT CHECKING
           if nargin < 2
               fileLabel = [];
           end 
           
           % CONFIRM PATH CONVENTION
           if ~strcmp(outputPath(end),'\')
               outputPath = strcat(outputPath,'\');
           end
                
           % GENERATE THE MONTE-CARLO SUBDIRECTORY
           timeStamp = datestr(datetime('now'),'yyyy-mm-dd@HH-MM-SS');              % Record current time
           fileString = sprintf('[%s] sessiondata%s',timeStamp,fileLabel);  % Build filestring
           sessionPath = strcat(outputPath,fileString);
           
           % CREATE OUTPUT DIRECTORY
           assert(mkdir(sessionPath) == 1,'Unable to create session directory.');
       end
       % BUILD THE RELATIVE MONTE-CARLO DIRECTORY
       function [monteCarloPath] = createMonteCarloDirectory(outputPath,fileLabel)
           % This function build the monte-carlo data output directory. The
           % naming of the file is of the form:
           % 'monteCarlo sessiondata 2018-09-28 Aobjects Bcycles
           
           % INPUT CHECKING
           if nargin < 2
               fileLabel = [];
           end 
           
           % CONFIRM PATH CONVENTION
           if ~strcmp(outputPath(end),'\')
               outputPath = strcat(outputPath,'\');
           end
                
           % GENERATE THE MONTE-CARLO SUBDIRECTORY
           timeStamp = datestr(datetime('now'),'yyyy-mm-dd@HH-MM-SS');              % Record current time
           fileString = sprintf('[%s] monteCarloInstance %s',timeStamp,fileLabel);   % Build filestring
           monteCarloPath = strcat(outputPath,fileString);
           
           % CREATE OUTPUT DIRECTORY
           assert(mkdir(monteCarloPath) == 1,'Unable to create Monte-Carlo output directory.');
       end
       
       % ///////////////// MONTE-CARLO SETUP UTILITIES ////////////////////
       % CONFIGURE THE PARALLEL/MULTITHREADED WORKERS
       function [poolObject]    = threadpoolInterface(conf)
           % This function is designed to generate a threadpool for the forth comming
           % application if one does not already exist. This is only necessary if the
           % aagent loops are sufficiently complicated (i.e. dt > 1s)
           % OUTPUTS:
           % poolObj - The pool object
           
           switch upper(char(conf))
               case 'CONFIRM'
                   fprintf('[MONTE-CARLO]\tConfirming parallel pool...\n');
                   % If exists, dont recreate.
                   poolObject = gcp('nocreate'); % Get current pool status
                   if isempty(poolObject)
                       % Open new pool for the monte-carlo simulations
                       simCluster = parcluster('local');
                       % CREATE THE PARALLEL POOL OBJECT
                       poolObject = parpool(simCluster,'IdleTimeout',30,'SpmdEnabled',true); % Generate thread pool, 2hour timeout         
                   end
                case 'RENEW'
                   fprintf('[MONTE-CARLO]\tRenewing parallel pool...\n');
                   % TERMINATE THE POOL OBJECT
                   delete(gcp('nocreate'));
                   
                   % Open new pool for the monte-carlo simulations
                   simCluster = parcluster('local');
                   % CREATE THE PARALLEL POOL OBJECT
                   poolObject = parpool(simCluster,'IdleTimeout',30,'SpmdEnabled',true); % Generate thread pool, 2hour timeout         
                   fprintf('[MONTE-CARLO]\t... Done.\n');                 
               case 'KILL'
                   fprintf('[MONTE-CARLO]\tDeleting current active pool...\n');
                   delete(gcp('nocreate'));
                   fprintf('[MONTE-CARLO]\t... Done.\n');
                   return
               otherwise
                   warning('Unrecognised threadpool interface command.');
                   return
           end
           fprintf('[MONTE-CARLO]\t... parallel pool ready.\n');
       end
       % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
       function [parsedConfig]  = configurationParser(defaultConfig,inputParameters)
           % This function is designed to parse a generic set of user
           % inputs and allow them to be compared to a default input
           % structure. This should be called in the class constructor
           
           % IF THE INPUT PARAMETERS ARE PROVIDED AS A STRUCTURE
           if isstruct(inputParameters)
               fieldLabels = fieldnames(inputParameters);
               parameterVector = cell(2*numel(fieldLabels),1);
               for i = 1:numel(fieldLabels)
                   parameterVector{2*i - 1} = fieldLabels{i};
                   parameterVector{2*i} = inputParameters.(fieldLabels{i});
               end
               inputParameters = parameterVector;
           end
           
           assert(mod(numel(inputParameters),2) == 0,'There are an uneven number of parameter-value pairs.');
           
           % ASSIGN CONFIG TEMPLATE
           parsedConfig = defaultConfig;
           pairNum = numel(inputParameters)/2;
           for n = 1:pairNum
               % PULL PARAMETER/VALUE PAIRS
               parameterLabel = inputParameters{2*n - 1};
               parameterValue = inputParameters{2*n};
               
               if isstruct(parsedConfig)
                   % IF THE 'TIME' SUBSTRUCTURE HAS THAT PROPERTY
                   if isfield(parsedConfig,parameterLabel)
                       parsedConfig.(parameterLabel) = parameterValue;
                   end
               else
                   % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                   if isprop(parsedConfig,parameterLabel)
                       parsedConfig.(parameterLabel) = parameterValue;   % Make a substitution
                   end
               end
           end
       end
       % GET THE DEFAULT MONTE-CARLO SETTINGS
       function [defaultConfig] = getDefaultMCconfig()
           % This function generates the default configuration for the
           % Monte-Carlo instance and the META parameters for all
           % subsequent OpenMAS instances.
           
           [~, userdir] = system('echo %USERPROFILE%');                    % Get desktop path
           
           % DEFAULT (GENERIC) MONTE-CARLO SETTINGS
           defaultConfig = struct();
           defaultConfig.outputPath = strcat(userdir,'\desktop\US18_data');  
           % DEFAULT MONTE-CARLO PARAMETERS 
           defaultConfig.objects = {};                                     % The object set
           defaultConfig.cycles = 1;                                       % Default number of cycles
           defaultConfig.phase = 'MONTE-CARLO';
           defaultConfig.threadPool = logical(false);                      % Complete in a asynchronous way
           defaultConfig.shutDownOnComplete = logical(false);              % Shutdown on completion
           defaultConfig.outputPath = pwd;
       end
       % GET THE DEFAULT OMAS-INSTANCE SETTINGS
       function [defaultConfig] = getDefaultOMASConfig()
           % This function creates a default configuration for any OpenMAS
           % configuration called by the OMAS_monteCarlo function.
           % INFO:
           % "agents/objects" - Not set initially.
           
           % DEFAULT (GENERIC) OpenMAS SETTINGS
           defaultConfig = struct();
           % HIGH-LEVEL I/O OPERATIONS
           defaultConfig.monteCarloMode = logical(true);
           defaultConfig.threadPool = logical(false);
           defaultConfig.verbosity = 0;                                    % Append simulation parameters
           defaultConfig.gui = logical(false);
           defaultConfig.figures = {'none'};
           % EVENT CONDITIONS
           defaultConfig.warningDistance = 10;
           defaultConfig.conditionTolerance = 1E-3;
           defaultConfig.visabilityModifier = 1;
           defaultConfig.positionSigma = 0;
           % TIMING CONFIGURATIONS
           defaultConfig.startTime = 0;
           defaultConfig.endTime = 10;
           defaultConfig.duration = defaultConfig.endTime - defaultConfig.startTime;
           defaultConfig.dt = 0.1;
       end
   end
end

