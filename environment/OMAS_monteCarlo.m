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
        objects = {};                           % The object set
        cycles = 1;                             % Default number of cycles
        sessions;                               % Parameter representing the num
        sessionLabels = {};                     % Parameter providing reference to each session
        isParallel = false;                     % Complete in a asynchronous way
        shutDownOnComplete = false;             % Shutdown on completion
        directory = pwd;
    end
    properties (Access = private)
        phase = 'MC-SETUP';
        OMASFiles = {'META','DATA','EVENTS','OBJECTS'};
        maxObjectNumber = 600;                  % The maximum number of objects defining studies
        sessionMETA;
        sessionID = 1;
    end    
    % ////////////////////////////// MAIN /////////////////////////////////
    methods
        % Constructor
        function [obj] = OMAS_monteCarlo(varargin)
            % ///////////////// Default settings //////////////////////////
            % Get the system paths
            pathString = mfilename('fullpath');
            pathString = pathString(1:(strfind(pathString,'environment')-1));
            addpath(OMAS_system.GetOSPathString([pathString,'environment\common'])); % Add system paths
            % Get the default output path
            [~, userdir] = system('echo %USERPROFILE%');                    % Get desktop path
            obj.directory = OMAS_system.GetOSPathString([userdir,'\desktop\US18_data']);
            
            % Get the default OpenMAS settings
            [obj.OMAS_settings] = obj.DefineDefaultOMASConfig();
                     
            % //////////////// Check for user overrides ///////////////////
            % Parse the overrides for the MC settings
            [obj] = GetParameterOverrides_recursive(obj,varargin);
            % /////////////////////////////////////////////////////////////
            
            % CONFIRM PROPOSED OBJECTS ARE VALID
            assert(iscell(obj.objects),'Please provide a object array of the form: [numberOfAlgorithms x numberOfObjects]');
            assert(numel(obj.objects) <= obj.GetMaxObjects(),sprintf('The number objects (%d) in the proposed study is to high.',numel(obj.objects)))
            
            % ////////////////// Post-initialisation //////////////////////            
            obj.sessions  = size(obj.objects,1);                           % Record the number of sessions
            % Generate session labels
            for i = 1:obj.sessions
                obj.sessionLabels{i} = class(obj.objects{i,1});
            end
            obj = obj.DefineSessionArray();                                % Build the descriptor array
            obj.directory = OMAS_system.GetOSPathString(obj.directory);    % Attempt to OS-normalise the requested paths
            
            % //////////////////////////////////////////////////////////////
            fprintf('\n[%s]\tIntialising Monte-Carlo instance.\n',obj.GetPhase);
            % ENABLE MATLAB PAUSING
            pause on
        end
    end
    methods
        % EVALUATE THE PROPOSED CYCLE
        function [obj] = EvaluateAllCycles(obj)
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
            obj.directory = obj.DefineDirectory(obj.directory,'monte_carlo_data'); % Define the path to the figure output
            obj.directory = OMAS_system.GetOSPathString([obj.directory,'\']);      % Change the directory    
            % IF REQUIRED, ENSURE A POOL IS PRESENT
            if obj.isParallel
                obj.ThreadpoolInterface('CONFIRM');
            end
            
            % ////////////// POPULATE THE SESSION DIRECTORIES //////////////
            obj = obj.SetPhase('CYCLING');
            for i = 1:obj.sessions
                % Each session represents a unique object set that will be
                % simulated over 'n' cycles.
                
                % //////////// DEFINE THE SESSION PARAMETERS ///////////////
                % Prepare the object vector
                sessionObjects = obj.objects(i,:);
                sessionObjects = sessionObjects(cellfun('isempty',sessionObjects) == 0);
                % Define session output directory
                sessionLabel = sprintf('-n(%d)-%s',numel(sessionObjects),obj.sessionLabels{i});
                sessionPath  = obj.DefineDirectory(obj.directory,['session_data',sessionLabel]);
                
                % ///////////////// COMPUTE THE CYCLES /////////////////////
                if obj.isParallel
                    obj = obj.ParallelCycles(sessionPath,sessionObjects);               % Compute session in parallel
                else
                    obj = obj.SequentialCycles(sessionPath,sessionObjects);             % Compute sessions sequentially
                end
                % //////////////////////////////////////////////////////////
                
                % Get the session descriptor
                [s] = obj.GetSessionMETA(i);
                
                % Print the session summary
                fprintf('\n[%s]\tSession complete; %d cycles successful (%d failed).\n',...
                    obj.GetPhase,sum(s.cycleLogicals),(obj.cycles - sum(s.cycleLogicals)));
                
                % Iterate the session
                obj = obj.SetSessionID(obj.GetSessionID + 1); 
                
                % Clean up
                clearvars -except obj highestLevelPath
            end
            
            % /////////////////// Hand over to analysis ////////////////////
            obj = obj.SetPhase('ANALYSIS');         % Define next phase
            obj = obj.SetSessionID(1);              % Reset the session number       
            
            % //////////////////// DATA IMPORTATION ////////////////////////
            % The number of data points imported can potentially be very
            % large. For this reason, the parallel pool is retained.
            % The data that will be imported.
            
            % Get the complete session META structure
            METAseries = obj.GetSessionMETA();
            if obj.isParallel
                parfor n = 1:obj.sessions
                    % COMPUTE SESSION BASED STATISTICS
                    METAseries(n) = obj.SessionAnalyitics(METAseries(n));
                end
            else
                for n = 1:obj.sessions
                    % COMPUTE SESSION BASED STATISTICS
                    METAseries(n) = obj.SessionAnalyitics(METAseries(n));
                end
            end
            % Override the session META with the new values
            obj = obj.SetSessionMETA(METAseries);
            
            % /////// CALCULATE THE MONTE-CARLO INSTANCE STATISTICS ////////
            % Using the data collected from each of the monte-carlo sessions
            % we can calculate and compare statistics from all sessions.
            obj = obj.SetPhase('OUTPUT');
                        
            % GET THE OpenMAS FIGURE PROPERTIES STRUCTURE (for continuity)
            [figureProperties] = OMAS_figureProperties();
            
            % GENERATE THE MEAN COMPUTATIONAL TIMES
            obj.GetMeanComputationTimeSeries(figureProperties,obj.directory,obj.GetSessionMETA());
            % GENERATE THE FIGURES COMPARING MEAN COMPUTATIONAL TIMES TO AGENT POPULATION
            obj.GetMeanComputationTimesVsPopulation(figureProperties,obj.directory,obj.GetSessionMETA());
            % GENERATE THE FIGURES COMPARING MEAN NUMBER OF COLLISIONS TO AGENT POPULATION
            obj.GetMeanCollisionsVsPopulation(figureProperties,obj.directory,obj.GetSessionMETA());
            % //////////////////////////////////////////////////////////////
            
            % KILL THE POOL WHEN ALL SESSIONS/STUDIES ARE COMPLETE
            if obj.isParallel
                obj.ThreadpoolInterface('KILL');
            end
        end
        % COMPUTE ANALYSIS OVER A SESSIONS - CYCLE SET
        function [sessionMETA] = SessionAnalyitics(obj,sessionMETA)
            % This function executes a series of basic analyses based on the
            % populated cycle directories. This is essentially the mean data
            % from the cycles described by the session conditions.
            
            % Sanity check
            if sum(sessionMETA.cycleLogicals) == 0
                fprintf('[%s]\tNo valid cycles detected, aborting session analysis.\n',obj.GetPhase());
                return
            end
            
            % Notes:
            % - The order of the files in the directory is changed due to
            %   windows. The cycle logicals do not correspond to the
            %   directories.
            
            % Get the cycle directories
            cycleDirectories = obj.GetPathSubDirectories(sessionMETA.path);
            totalSamples = numel(cycleDirectories);
            validSamples = sum(sessionMETA.cycleLogicals);
            sessionMETA.errors = totalSamples - validSamples;              % Errors
            MEANS = sessionMETA.MEANS;                                     % Get a local copy of the MEANS structure      
                        
            % CREATE THE SUMMARY CSV
            fileName = OMAS_system.GetOSPathString([sessionMETA.path,'\','summaryData.csv']);
            obj.PushDataToCSV(fileName,1,{'collisions','waypoints','mean_dt','min_dt','max_dt','minSeparation','maxSeparation'});
            
            % GET THE DATA FROM THE FIRST CYCLE
            [isSuccessful,initialSample] = obj.ImportPathData(obj.GetOMASFiles,[sessionMETA.path,'\',cycleDirectories{1}]);
            if any(isSuccessful == false)
                warning('Unable to read first cycle data, aborting');
                return
            end
            
            % GET THE NUMBER OF AGENTS IN THE SIMULATION
            % This is only valid if all cycles are computed with the same
            % OMAS_settings, which they should be.
            sessionMETA.agents     = initialSample.META.totalAgents;
            sessionMETA.waypoints  = initialSample.META.totalWaypoints;            
            sessionMETA.timeVector = initialSample.META.TIME.timeVector;
            clearvars initialSample
            
            % /////////////// BEGIN RUNNING THROUGH CYCLES /////////////////
            for i = 1:validSamples
                % The path to cycles in the 'valid' subset
                cyclePath = OMAS_system.GetOSPathString([sessionMETA.path,'\',cycleDirectories{i}]);
                
                % IMPORT THE CYCLE SAMPLE
                [isSuccessful,cycleSample] = obj.ImportPathData(obj.GetOMASFiles,cyclePath);
                if any(isSuccessful == false)
                    fprintf('[%s]\tSkipping failed cycle in: \n[%s]\tSession: "%s" \n[%s]\tCycle: "%s"\n',...
                        obj.GetPhase(),obj.GetPhase(),sessionMETA.path,obj.GetPhase(),cycleDirectories{i}); 
                    continue
                end
                
                % //////////////// COMPUTE CYCLE DATA /////////////////////
                % COMPUTE THE AGENT-TIMING STATISTICS
                [mean_dt,min_dt,max_dt,dt_timeSeries,~] = obj.cycle_timingAnalytics(...
                    cycleSample.META,...
                    cycleSample.OBJECTS);
                % COMPUTE THE AGENT-TRAJECTORY STATISTICS
                [minSeparation,maxSeparation] = obj.cycle_trajectoryAnalytics(...
                    cycleSample.META,...
                    cycleSample.DATA);
                % COMPUTE THE CYCLE - EVENT BASED STATISTICS
                [collisions,waypoints] = obj.cycle_eventAnalytics(...
                    cycleSample.DATA);
                
                % Push the data to a .CSV (for reference)
                obj.PushDataToCSV(fileName,0,[collisions,waypoints,mean_dt,min_dt,max_dt,minSeparation,maxSeparation])
                
                % /////////////// UPDATE MEAN VALUES //////////////////////
                % Mean events
                MEANS.collisions = MEANS.collisions + collisions;
                MEANS.waypoints  = MEANS.waypoints + waypoints;
                
                % Separation (min/max) 
                MEANS.minSeparation = MEANS.minSeparation + minSeparation;
                MEANS.maxSeparation = MEANS.maxSeparation + maxSeparation;
                
                % Computation timing coefficients
                MEANS.loopTime    = MEANS.loopTime + mean_dt;
                MEANS.minLoopTime = MEANS.minLoopTime + min_dt;
                MEANS.maxLoopTime = MEANS.maxLoopTime + max_dt;
                
                % Define the dt time-series
                if i == 1
                    MEANS.dt_timeSeries = dt_timeSeries;    
                else
                    MEANS.dt_timeSeries = MEANS.dt_timeSeries + dt_timeSeries; 
                end
            end
            
            % Calculate the means
            fn = fieldnames(MEANS);
            for i = 1:length(fn)
                MEANS.(fn{i}) = MEANS.(fn{i})./validSamples;
            end
            
            % Reassign for structure
            sessionMETA.MEANS = MEANS;
            
            % //////////// PUSH SESSION MEAN DATA TO .mat FILE /////////////
            % The cycle path
            save(OMAS_system.GetOSPathString([sessionMETA.path,'\MEANS.mat']), '-struct', 'sessionMETA');
            
            % Clear trash
            clearvars -except sessionMETA;
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
            mean_dt = 0; max_dt  = 0; min_dt  = 0;
            dt_timeSeries = zeros(size(timeVector));
            
            % CALCULATE THE AGENT MEANS
            for agentNum = 1:agentNumber
                % DATA LOCAL TO A GIVEN AGENT OF A GIVEN CYCLE
                agentDATA = agentIndex{agentNum}.DATA;
                
                % Sanity check
                if ~isfield(agentDATA,'indicator') || ~isfield(agentDATA,'dt')
                    continue
                end
                
                % DETERMINE THE ALGORITHM COMPUTATION TIME TIME-SERIES
                valid_algorithm_dt = agentDATA.dt(logical(agentDATA.indicator));                   % Get the times where the computations were ran
                % GET THE (AGENT) ALGORITHM COMPUTATION-TIME DATA
                agent_mean_dt = sum(valid_algorithm_dt)/numel(valid_algorithm_dt);
                agent_max_dt  = max(valid_algorithm_dt);
                agent_min_dt  = min(valid_algorithm_dt);
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
            
            % Sanity check
            assert(isstruct(cycle_DATA),'The sample does not have a DATA sub-structure.');
            
            % GET THE NUMBER OF CYCLE COLLISIONS
            collisions = double(cycle_DATA.collisions);
            % GET THE NUMBER OF CYCLE WAY-POINTS ACHIEVED
            waypoints  = double(cycle_DATA.waypointsAchieved);
        end
    end
    %% ///////////////////// DATA HANDLING UTILITIES //////////////////////
    methods (Static)
        % ////////////////////// FIGURE GENERATION /////////////////////////
        % GET THE COMPARISON OF COLLISIONS TO POPULATION
        function [figureHandle] = GetMeanCollisionsVsPopulation(figureProperties,figurePath,meanData)
            % This function computes the mean number of collisions against
            % on population.
            
            % CONFIGURE THE PLOT ATTRIBUTES
            figurePath = strcat(figurePath,'mean_collisions_population');
            figureHandle = figure('Name','OpenMAS Collisions vs Population');
            set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
            set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
            ax = axes(figureHandle);
            
            % REORDER THE DATA BY NUMBER OF AGENTS
            [~,orderIndices] = sort([meanData.agents]);
            meanData = meanData(orderIndices);          % Reorder mean-data
            
            plottedLabels = {};
            for i = 1:numel(meanData)
                % No session data
                if sum([meanData(i).cycleLogicals]) == 0
                    continue
                end
                % Check if the data series has already been plotted
                if any(ismember(plottedLabels,meanData(i).label))
                    continue
                else
                    plottedLabels = horzcat(plottedLabels,meanData(i).label);
                end

                % Get the mean data where 'label' is common
                plotData = meanData(ismember({meanData.label},meanData(i).label));
                collisionSeries  = zeros(size(plotData));
                for j = 1:length(plotData)
                   collisionSeries(j) = plotData(j).MEANS.collisions;    % The mean minimum time-series
                end
                populationSeries = [plotData(:).agents]';                % The population series
                
                % DEFINE ERROR BAR PLOT
                hold on;
                plot(ax,populationSeries,collisionSeries,...
                    'LineWidth',figureProperties.lineWidth,...
                    'Color',rand(3,1),...
                    'DisplayName',['$',strrep(plotData(1).label,'_','-'),'$'],...
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
                'Agent population',...
                'interpreter',figureProperties.interpreter,...
                'fontweight',figureProperties.fontWeight,...
                'fontSize',figureProperties.axisFontSize,...
                'FontSmoothing','on');
%             xlim([0,meanData.agents]);
            % Y-Label
            ylabel(ax,...
                'Mean collisions',...
                'interpreter',figureProperties.interpreter,...
                'fontweight',figureProperties.fontWeight,...
                'fontSize',figureProperties.axisFontSize,...
                'FontSmoothing','on');
%             ylim([0,meanData.agents]);
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
                'String',sprintf('Cycles: %d',sum(meanData(1).cycleLogicals)),...
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
        function [figureHandle] = GetMeanComputationTimesVsPopulation(figureProperties,figurePath,meanData)
            % This function computes the mean computation times of each
            % session against the number of agents in its population.
            
            % CONFIGURE THE PLOT ATTRIBUTES
            figurePath = strcat(figurePath,'mean_computation_time_population');
            figureHandle = figure('Name','OpenMAS Computation Time vs Population');
            set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
            set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
            ax = axes(figureHandle);
            
            % REORDER THE DATA BY NUMBER OF AGENTS
            [~,orderIndices] = sort([meanData.agents]);
            meanData = meanData(orderIndices);          % Reorder mean-data
            
            plottedLabels = {};
            for i = 1:numel(meanData)
                % No session data
                if sum([meanData(i).cycleLogicals]) == 0
                    continue
                end
                % Check if the data series has already been plotted
                if any(ismember(plottedLabels,meanData(i).label))
                    continue
                else
                    plottedLabels = horzcat(plottedLabels,meanData(i).label);
                end
                
                % /////////////////// Generate data ///////////////////////
                % Get the mean data where 'label' is common
                plotData = meanData(ismember({meanData.label},meanData(i).label));
                minSeries  = zeros(size(plotData));
                meanSeries = zeros(size(plotData));
                maxSeries  = zeros(size(plotData));
                for j = 1:length(plotData)
                   minSeries(j)  = plotData(j).MEANS.minLoopTime;    % The mean minimum time-series
                   meanSeries(j) = plotData(j).MEANS.loopTime;       % The mean time-series
                   maxSeries(j)  = plotData(j).MEANS.maxLoopTime;    % The mean maximum times-series
                end
                populationSeries = [plotData(:).agents]';            % The population series
                
                % DEFINE ERROR BAR PLOT
                hold on;
                errorbar(ax,populationSeries,meanSeries,minSeries,maxSeries,...
                    'LineWidth',figureProperties.lineWidth,...
                    'Color',rand(3,1),...
                    'DisplayName',['$',strrep(plotData(1).label,'_','-'),'$'],...
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
                'Agent population',...
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
%             ylim([0,meanData.MEANS.maxLoopTime]);
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
                'String',sprintf('Cycles: %d',sum(meanData(1).cycleLogicals)),...
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
        function [figureHandle] = GetMeanComputationTimeSeries(figureProperties,figurePath,meanData)
            % This function computes the mean computation time-series of all
            % agents and compares them with each session.
            
            % ASSUMPTIONS:
            % meanData - The mean times are given in seconds.
            
            % CONFIGURE THE PLOT ATTRIBUTES
            figurePath = strcat(figurePath,'mean_computation_times');
            figureHandle = figure('Name','OpenMAS Mean Computational Timeseries');
            set(figureHandle,'Position',figureProperties.windowSettings);         % [x y width height]
            set(figureHandle,'Color',figureProperties.figureColor);               % Background colour
            ax = axes(figureHandle);
            
            for sessionNum = 1:length(meanData)
                % No session data
                if sum([meanData.cycleLogicals]) == 0
                    continue
                end
                hold on;
                % NAME BASED ON SESSION
                displayName = ['n=',num2str(meanData(sessionNum).agents),',$ $',meanData(sessionNum).label];
                % PLOT THE MEAN TRACES FOR EACH SESSION
                stairs(ax,...
                    meanData(sessionNum).timeVector',...
                    meanData(sessionNum).MEANS.dt_timeSeries',...
                    'Color',rand(3,1),...
                    'LineWidth',figureProperties.lineWidth,...
                    'DisplayName',['$',strrep(displayName,'_','-'),'$']);

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
                'String',sprintf('Cycles: %d',sum(meanData(1).cycleLogicals)),...
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
        
        % IMPORT OMAS FILES FROM A SPECIFIC SESSION DIRECTORY
        function [successVector,cycleSample] = ImportPathData(OMASFileNames,cyclePath)
            % This function imports the specified ".mat" files from the
            % specified path. This is typically the output path of a given
            % OpenMAS cycle.
                        
            % Sanity check
            assert(ischar(cyclePath),'Please provide a valid path to the session directory');
            assert(iscell(OMASFileNames),'Expecting a cell array of file names');
            
            % Constants
            cycleSample = cell2struct(cell(size(OMASFileNames)),OMASFileNames,2);
            successVector = true(length(OMASFileNames),1);
  
            % IMPORT THE SESSION .MAT FILES
            for fileName = 1:numel(OMASFileNames)
                % PATH TO OMAS VARIABLES
                filePath = strcat(cyclePath,'\',OMASFileNames{fileName},'.mat');
                % Parse the path into a string friendly to the host OS
                filePath = OMAS_system.GetOSPathString(filePath);
                % Default return
                cycleSample.(OMASFileNames{fileName}) = [];
                try
                    load(filePath,'-mat');
                    % IMPORT FILE AS VARIABLE
                    cycleSample.(OMASFileNames{fileName}) = eval(OMASFileNames{fileName});
                    % Log string
%                     logString = horzcat(logString,[OMASFileNames{fileName},' imported correctly. ']);
                catch importError
%                     logString = horzcat(logString,[OMASFileNames{fileName},': ',importError.message,' ']);
                    successVector(fileName) = false;
                end
            end  
        end
        % PUSH VALUES TO CSV FILE
        function PushDataToCSV(filename,isHeaderLogical,rowData)
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
        % PUSH MESSAGE TO LOG
        function PushMessageToLog(sessionPath,message)
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
    end
    
    % ////////////////////////// CYCLE FUNCTIONS //////////////////////////
    methods
        % RUN A DEFINED NUMBER OF SEQUENCIAL CYCLES
        function [obj] = SequentialCycles(obj,sessionPath,objectIndex)
            % This function computes a monte-carlo sequence of OMAS
            % instances.
            
            % OPEN LOG (Only one handle is needed)
            logPath  = [sessionPath,'\log.txt'];
            simPhase = obj.GetPhase();
            sessionData = obj.GetSessionMETA(obj.GetSessionID);
            sessionData.path = sessionPath;
            isSuccessfulArray = sessionData.cycleLogicals;
            
            % /////////////// BEGIN THE CYCLE ITERATIONS ///////////////////
            fprintf('\n[%s]\tEXECUTING %d SEQUENTIAL MONTE-CARLO CYCLES...\n\n',simPhase,obj.cycles);
            for cycle = 1:obj.cycles
                fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',simPhase,cycle);
                [isSuccessfulArray(cycle),cycleMessage] = obj.SingleCycle(sessionPath,obj.OMAS_settings,objectIndex);
                fprintf('\n[%s]\tCYCLE %d COMPLETE\n\n',simPhase,cycle);
                % PUSH CYCLE SUMMARY TO LOG FILE
                obj.PushMessageToLog(logPath,sprintf('(Cycle:%d)\t%s',cycle,cycleMessage));
            end
            fprintf('[%s]\t%d SEQUENTIAL CYCLES COMPLETE.\n',simPhase,obj.cycles);
            % //////////////////////////////////////////////////////////////
            
            % Post-processing
            sessionData.cycleLogicals = isSuccessfulArray;                 % Record the cycle successes
            obj = obj.SetSessionMETA(sessionData,obj.GetSessionID);        % Record the session statistics
            clearvars -except obj                                          % Clear trash
        end
        % RUN A DEFINED NUMBER OF PARALLEL CYCLES
        function [obj] = ParallelCycles(obj,sessionPath,objectIndex)
            % Compute the proposed number of cycles in a parallel manner.
            % The resulting data generated by each simulation instance is
            % output to each session directory.
            
            fprintf('[%s]\tConfirming parallel worker pool...\n',obj.GetPhase());
            
            % Preparation
            obj.ThreadpoolInterface('CONFIRM');                         % Confirm the pool
            logPath  = [sessionPath,'\log.txt'];
            simPhase = obj.GetPhase();
            sessionData = obj.GetSessionMETA(obj.GetSessionID);
            sessionData.path = sessionPath;
            isSuccessfulArray = sessionData.cycleLogicals;
            
            % /////////////// BEGIN THE CYCLE ITERATIONS //////////////////
            fprintf('\n[%s]\tEXECUTING %s PARALLEL MONTE-CARLO CYCLES...\n\n',simPhase,num2str(obj.cycles));
            parfor cycle = 1:obj.cycles
                fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',simPhase,cycle);
                [isSuccessfulArray(cycle),cycleMessage] = obj.SingleCycle(sessionPath,obj.OMAS_settings,objectIndex);
                fprintf('\n[%s]\tCYCLE %d COMPLETE.\n\n',simPhase,cycle);
                % PUSH CYCLE SUMMARY TO LOG FILE
                obj.PushMessageToLog(logPath,sprintf('(Cycle:%d)\t%s',cycle,cycleMessage));
            end
            fprintf('[%s]\t%d PARALLEL CYCLES COMPLETE.\n',simPhase,obj.cycles);
            % /////////////////////////////////////////////////////////////
            
            % Post-processing
            sessionData.cycleLogicals = isSuccessfulArray;                 % Record the cycle successes
            obj = obj.SetSessionMETA(sessionData,obj.GetSessionID);        % Record the session statistics
            clearvars -except obj                                          % Clear trash
        end
    end
    methods (Static)
        % RUN ONE PARAMETERISED SIMULATION CYCLE
        function [isSuccessful,cycleMessage] = SingleCycle(sessionPath,SETTINGS,objectIndex)
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
                    objectIndex{objectNum}.SetGLOBAL('position',...
                        + objectIndex{objectNum}.GetGLOBAL('position')...
                        + [SETTINGS.positionSigma*randn(2,1);0]);          	% 2D perturbation;
                    % PERTURB THE GLOBAL VELOCITIES
                    objectIndex{objectNum}.SetGLOBAL('velocity',...
                        + objectIndex{objectNum}.GetGLOBAL('velocity')...
                        + [SETTINGS.velocitySigma*randn(2,1);0]);          	% 2D perturbation;
                catch perturbationError
                    warning('Unable to perturb global position of object %s',num2str(objectNum));
                    rethrow(perturbationError);
                end
            end
            clearvars -except sessionPath SETTINGS objectIndex
            % //////////////////////////////////////////////////////////////
            
            % ////// CREATE OpenMAS INSTANCE WITH THE DEFINED SETTINGS /////
            try
                OMAS_initialise(...
                    'idleTimeOut',SETTINGS.idleTimeOut,...
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
                isSuccessful = true;
            catch cycleError
                % FAILURE MESSAGE
                warning('\n...cycle error.\n');
                cycleMessage = cycleError.message;
                disp(getReport(cycleError, 'extended', 'hyperlinks', 'on' ));
                isSuccessful = false;
            end
            clearvars -except isSuccessful cycleMessage;
            % //////////////////////////////////////////////////////////////
            
            % PAUSE FOR STABILITY
            pause(0.1); % 0.1s delay
        end
    end
    
    % /////////////////////////// SETUP METHODS ///////////////////////////
    methods
        % Define session descriptor array
        function [obj] = DefineSessionArray(obj)
            % Create the descriptor structure
            s = struct();
            s.label = [];
            s.path = [];
            s.cycleLogicals = [];
            s.errors = 0;
            
            % The number of agents and waypoints
            s.agents = 0;
            s.waypoints = 0;
            s.timeVector = NaN(1,2);
            
            % Mean-data substructure
            s.MEANS = struct();
            s.MEANS.collisions = 0;
            s.MEANS.waypoints = 0;
            s.MEANS.minSeparation = 0;
            s.MEANS.maxSeparation = 0;
            s.MEANS.loopTime = 0;
            s.MEANS.minLoopTime = 0;
            s.MEANS.maxLoopTime = 0;
            s.MEANS.dt_timeSeries = NaN(1,2);
            
            % Concatinate
            s = repmat(s,[obj.sessions, 1]);
            for i = 1:obj.sessions
                s(i).label = obj.sessionLabels{i};
                s(i).cycleLogicals = false(obj.cycles,1);
            end            
            % Assign the array to the object
            obj = obj.SetSessionMETA(s);
        end 
    end
     methods (Static)
        % GET LIST OF SUB-DIRECTORIES
        function [subDirectories] = GetPathSubDirectories(pathString)
            % INPUT CHECK
            assert(ischar(pathString),'Please provide a valid path string.');
            % FETCH SESSION FOLDER SET
            rawDirectories = dir(pathString);
            % OMIT KNOWN BAD PATHS
            subDirectories = rawDirectories(~ismember({rawDirectories.name},{'.','..'}));
            % OMIT NON-FOLDERS ( JUST THE LIST OF SESSION DIRECTORIES )
            subDirectories = subDirectories([subDirectories.isdir] == 1);
            % GET THE PATH STRINGS
            subDirectories = {subDirectories.name}';
        end
        % BUILD THE RELATIVE MONTE-CARLO DIRECTORY
        function [monteCarloPath] = DefineDirectory(outputPath,fileLabel)
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
            timeStamp = datestr(datetime('now'),'yyyy-mm-dd @ HH-MM-SS');              % Record current time
            fileString = sprintf('[%s] %s',timeStamp,fileLabel);   % Build filestring
            monteCarloPath = strcat(outputPath,fileString);
            
            % CREATE OUTPUT DIRECTORY
            assert(mkdir(monteCarloPath) == 1,'Unable to create directory.');
        end
        % CONFIGURE THE PARALLEL/MULTITHREADED WORKERS
        function [poolObject]     = ThreadpoolInterface(conf)
            % This function is designed to generate a threadpool for the forth comming
            % application if one does not already exist. This is only necessary if the
            % aagent loops are sufficiently complicated (i.e. dt > 1s)
            % OUTPUTS:
            % poolObj - The pool object
            
            switch upper(char(conf))
                case 'CONFIRM'
                    fprintf('Confirming parallel pool...\n');
                    % If exists, dont recreate.
                    poolObject = gcp('nocreate'); % Get current pool status
                    if isempty(poolObject)
                        % Open new pool for the monte-carlo simulations
                        simCluster = parcluster('local');
                        % CREATE THE PARALLEL POOL OBJECT
                        poolObject = parpool(simCluster,'IdleTimeout',30,'SpmdEnabled',true); % Generate thread pool, 2hour timeout
                    end
                case 'RENEW'
                    fprintf('Renewing parallel pool...\n');
                    % TERMINATE THE POOL OBJECT
                    delete(gcp('nocreate'));
                    
                    % Open new pool for the monte-carlo simulations
                    simCluster = parcluster('local');
                    % CREATE THE PARALLEL POOL OBJECT
                    poolObject = parpool(simCluster,'IdleTimeout',30,'SpmdEnabled',true); % Generate thread pool, 2hour timeout
                    fprintf('... Done.\n');
                case 'KILL'
                    fprintf('Deleting current active pool...\n');
                    delete(gcp('nocreate'));
                    fprintf('... Done.\n');
                    return
                otherwise
                    warning('Unrecognised threadpool interface command.');
                    return
            end
            fprintf('... parallel pool ready.\n');
        end
        % GET THE DEFAULT OMAS-INSTANCE SETTINGS
        function [defaultConfig]  = DefineDefaultOMASConfig()
            % This function creates a default configuration for any OpenMAS
            % configuration called by the OMAS_monteCarlo function.
            % INFO:
            % "agents/objects" - Not set initially.
            
            % DEFAULT (GENERIC) OpenMAS SETTINGS
            defaultConfig = struct();
            % HIGH-LEVEL I/O OPERATIONS
            defaultConfig.monteCarloMode = true;
            defaultConfig.threadPool = false;
            defaultConfig.verbosity = 0;                                    % Append simulation parameters
            defaultConfig.gui = false;
            defaultConfig.figures = {'none'};
            % EVENT CONDITIONS
            defaultConfig.warningDistance = 10;
            defaultConfig.conditionTolerance = 1E-3;
            defaultConfig.visabilityModifier = 1;
            % TIMING CONFIGURATIONS
            defaultConfig.startTime = 0;
            defaultConfig.endTime = 10;
            defaultConfig.duration = defaultConfig.endTime - defaultConfig.startTime;
            defaultConfig.dt = 0.25;
            defaultConfig.idleTimeOut = 0;
            % Perturbations
            defaultConfig.positionSigma = 0;
            defaultConfig.velocitySigma = 0;
        end
    end   
    % ////////////////////////// GET/SET METHODS //////////////////////////
    methods
        % Get the current session ID
        function [n]   = GetSessionID(obj)
            n = obj.sessionID;
        end
        % Set the current session ID
        function [obj] = SetSessionID(obj,n)
            obj.sessionID = n;
        end
        % Get the session structure
        function [s]   = GetSessionMETA(obj,i)
            if nargin == 1
                i = ':';
            end
            s = obj.sessionMETA(i);
        end
        % Set the session structure
        function [obj] = SetSessionMETA(obj,s,i)
            if nargin < 3 
                i = ':';
            end
            assert(isstruct(s),'Expecting a structure.');
            
            if isempty(obj.sessionMETA)
                obj.sessionMETA = s;
            else
                obj.sessionMETA(i) = s;
            end
        end
        % Get the OMAS file names
        function [fn]  = GetOMASFiles(obj)
            fn = obj.OMASFiles;
        end
        % Get the maximum number of objects
        function [n]   = GetMaxObjects(obj)
            n = obj.maxObjectNumber;
        end    
        % Get the MC phase
        function [p]   = GetPhase(obj)
            p = obj.phase;
        end
        % Set the MC phase
        function [obj] = SetPhase(obj,p)
            obj.phase = p; 
        end
    end    
end

