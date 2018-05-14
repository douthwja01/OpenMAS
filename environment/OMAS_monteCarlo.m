%% OpenMAS MONTE-CARLO SIMULATION TOOLBOX (OMAS_monteCarlo.m) %%%%%%%%%%%%%
% This function provides a suite of tools for monte-carlo type simulations
% of defined agent sets and scenarios.

% Author: James A. Douthwaite 24/03/18

classdef OMAS_monteCarlo
   properties
       % DEFAULT SIMULATION PARAMETERS
       SETTINGS;
       % DEFAULT MONTE-CARLO PARAMETERS
       outputPath = strcat(pwd,'\data');                                   % Highest level output path
       systemFile = 'studyData.mat';
       sessionPath;                                                        % The directory of the individual study
       sessionNumber = 1; 
       cycles = 1;                                                         % Default number of cycles
       phase = 'MONTE-CARLO';
       OMASFileNames = {'META','DATA','EVENTS'};
       threadPool = logical(false);
       shutDownOnComplete = logical(false);
   end
   % ///////////////////// MONTE-CARLO UTILITIES //////////////////////////
   methods
       % CONSTRUCTOR METHOD
       function [obj] = OMAS_monteCarlo(varargin)
           % INPUT HANDLING
           if length(varargin) == 1 && iscell(varargin)                    % Catch nested cell array inputs
               varargin = varargin{:};
           end
           fprintf('\n[%s]\tIntialising Monte-Carlo instance.\n',obj.phase);
           % PARSE MONTE-CARLO PARAMETERS
           [obj] = obj.configurationParser(obj,varargin); 
                      
           % GET THE SIMULATION SETTINGS
           defaultConfig = obj.getDefaultConfig();  
           obj.SETTINGS  = obj.configurationParser(defaultConfig,varargin); % Parse input config against the default config
           
           % INPUT HANDLIING
           assert(~isempty(obj.SETTINGS.objects),'Please specify an object configuration cell array of the form: [numberOfAlgorithms x numberOfObjects]');
           
           % CONFIGURE THE MONTE-CARLO DIRECTORIES
           obj.sessionPath = obj.outputPath;
           obj.systemFile = [obj.sessionPath,'\',obj.systemFile];          % System files
           
           % ENABLE MATLAB PAUSING
           pause on
       end
   end   
   % ////////////////////////// DATA ANALYSIS /////////////////////////////
   methods
       % HIGHEST LEVEL 'N' CYCLE ANALYSIS
       function [summaryData,obj] = evaluateCycles(obj)
           % This function computes a complete cycle set, imports the
           % resultant data and computes the available statistical analyses.
           
           % ASSUMPTIONS:
           % - The .SETTINGS.objects input takes the form:
           %   [numberOfStudies x numberOfObjects]. This therefore
           %   allows sequential testings of multiple object sets.           
           
           %            % GENERATE WAIT-BAR
           %            wbHandle = waitbar(0,'Processing...',...
           %                'Name','Processing Monte-Carlo Cycles',...
           %                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
           %                'Position',[0;50;40;20],...
           %                'WindowStyle','normal');
           %            setappdata(wbHandle,'canceling',0);
           
           % DETERMINE THE NATURE OF THE PROCEDURES
           [numStudies,~] = size(obj.SETTINGS.objects);
           for study = 1:numStudies
               % DEFINE THE CURRENT STUDY NUMBER
               obj.sessionNumber = study;
               
               % GET THE NUMBER OF STUDY OBJECTS
               objectIndex = obj.SETTINGS.objects(obj.sessionNumber,:);
               numObjects = numel(objectIndex(cellfun('isempty',objectIndex) == 0)); % Omit empty cells
               clear objectIndex
               
               % ALLOCATE NEW SESSION PATH FOR THE GIVEN STUDY
               [obj.sessionPath] = obj.getMonteCarloDirectory(obj.outputPath,numObjects,obj.cycles);
               
               %                % CHECK WAIT-BAR INTERACTION
               %                if getappdata(wbHandle,'canceling')
               %                    delete(wbHandle);
               %                    break
               %                else
               %                    % UPDATE THE WAIT-BAR PROGRESS
               %                    wbStatus = sprintf('Session progress (%d/%d)',study,numStudies);
               %                    waitbar((study/numStudies),wbHandle,wbStatus);
               %                end
               
               % ///////////////////// RUN CYCLES /////////////////////////
               % BEGIN RUNNING CYCLES
               if obj.threadPool
                   [obj,lastCycle,cycleMessage] = obj.parallelCycles();
               else
                   [obj,lastCycle,cycleMessage] = obj.sequentialCycles();
               end
               
               if lastCycle ~= obj.cycles
                   logString = sprintf('Final cycle(%d/%d), session %d failed with error: \n %s',lastCycle,obj.cycles,obj.sessionNumber,cycleMessage);
                   obj.writeLogFile(obj.sessionPath,logString)
               else
                   logString = sprintf('Final cycle(%d/%d), session %d complete successfully.',lastCycle,obj.cycles,obj.sessionNumber);
                   obj.writeLogFile(obj.sessionPath,logString);
               end

               clearvars -except obj
               
               fprintf('\n[%s]\tAnalysing output data...\n',obj.phase);
               
               % GET THE VALID DATA PATH LIST
               [validSamplePaths] = obj.getValidSessionPaths(obj.sessionPath);
               numSamples = numel(validSamplePaths);
               
               % COLLATE MEAN DATA OVER SAMPLES
               summaryData = struct('meanCollisions',0,...
                                    'meanWaypoints',0,....
                                    'minSeparation',0,...
                                    'maxSeparation',0,...
                                    'meanLoopTime',0,...
                                    'minLoopTime',0,...
                                    'maxLoopTime',0,...
                                    'meanLoopTimeSeries',[],...
                                    'meanTimeVector',[]);
                                
               % MOVE THROUGH THE SESSION SAMPLES                 
               for sampleNum = 1:numSamples
                   % EVALUATE EVENT DATA AT A GIVEN SESSION PATH
                   absSamplePath = strcat(obj.sessionPath,'\',validSamplePaths(sampleNum));
                   [sample,isSuccessful] = obj.importPathData(absSamplePath);
                   if ~isSuccessful
                       continue
                   end
                   
                   % ////// UPDATE MONTE-CARLO DATA WITH NEW SAMPLE ///////
                   % UPDATE MEAN EVENT STATS
                   [collisions,waypoints] = obj.getEventStatistics(sample);
                   summaryData.meanCollisions = summaryData.meanCollisions + collisions/numSamples;
                   summaryData.meanWaypoints  = summaryData.meanWaypoints + waypoints/numSamples;
                   % UPDATE MEAN SEPARATION STATS
                   [minSeparation,maxSeparation] = obj.getTrajectoryStatistics(sample);
                   summaryData.minSeparation = summaryData.minSeparation + minSeparation/numSamples;
                   summaryData.maxSeparation = summaryData.maxSeparation + maxSeparation/numSamples;
                   % UPDATE MEAN AGENT/ALGORITHM TIMING PARAMETERS
                   [mean_dt,min_dt,max_dt,dt_timeSeries,timeVector] = obj.getAgentTimingStatistics(sample);
                   summaryData.meanLoopTime = summaryData.meanLoopTime + mean_dt/numSamples;
                   summaryData.minLoopTime = summaryData.minLoopTime + min_dt/numSamples;
                   summaryData.maxLoopTime = summaryData.maxLoopTime + max_dt/numSamples;
                   summaryData.meanTimeVector = timeVector;
                   if isempty(summaryData.meanLoopTimeSeries)
                       % IF THIS IS THE FIRST CYCLE-SAMPLE
                       summaryData.meanLoopTimeSeries = dt_timeSeries/obj.cycles;
                   else
                       % IF THE CURRENT MEAN VECTOR IS LONGER THAN THE CYCLE-SAMPLE
                       summaryData.meanLoopTimeSeries = summaryData.meanLoopTimeSeries + dt_timeSeries/numSamples;
                   end
                   % PURGE SAMPLE
                   clear sample absSamplePath
               end
               
               % ////////// GENERATE THE KNOWN OUTPUT FIGURES /////////////
               % THE AGENT COMPUTATIONAL TIME TIME-SERIES
               obj.get_MC_timeSeries(summaryData);

               % PUSH TO DATA TO OUTPUT DIRECTORY
               save(strcat(obj.sessionPath,'\summaryData.mat'),'summaryData');
               filePath = strcat(obj.sessionPath,'\summaryData.csv');
               writetable(struct2table(summaryData,'AsArray',true),filePath);
               
               % PUSH SUMMARY LABEL 
               obj.writeSessionDescriptor();
               
%                % CHECK WAIT-BAR INTERACTION
%                if getappdata(wbHandle,'canceling')
%                    delete(wbHandle);
%                    break
%                else
%                    % UPDATE THE WAIT-BAR PROGRESS
%                    wbStatus = sprintf('Session progress (%s/%s)',num2str(study),num2str(numStudies));
%                    waitbar(study/numStudies,wbHandle,wbStatus);
%                end
               
               % CLOSE FIGURES
               close all;
           end
           
           fprintf('\n[%s]\t..Exiting.\n',obj.phase);
           % SHUTDOWN IF REQUESTED
           if obj.shutDownOnComplete
               fprintf('\n[%s]\t...shutting down PC.\n',obj.phase);
               system('shutdown -s');
           end
           clearvars -except summaryData obj
       end
       % IMPORT A SESSION FROM PATH
       function [sessionSample,isSuccessful] = importPathData(obj,sessionPath)
           % This function is designed to import the OpenMAS session data
           % structures from a provided path as a more lean way of
           % importing data.          
           
           sessionSample = [];
           isSuccessful = 0;
           
           % INPUT HANDLING
           if iscell(sessionPath)
               sessionPath = char(sessionPath);
           elseif ~ischar(sessionPath)
               warning('Please provide a valid path to the session directory');
               return
           end

           % IMPORT THE SESSION .MAT FILES
           sessionSample = cell2struct(cell(size(obj.OMASFileNames)),obj.OMASFileNames,2);
           % FOR EACH KNOWN DATA OBJECT
           for OMASfile = 1:numel(obj.OMASFileNames)
               % PATH TO OMAS VARIABLES
               filePath = strcat(sessionPath,'\',obj.OMASFileNames{OMASfile},'.mat');
               try
                   load(filePath,'-mat');
                   % IMPORT FILE AS VARIABLE
                   sessionSample.(obj.OMASFileNames{OMASfile}) = eval(obj.OMASFileNames{OMASfile});
               catch
                   warning('Unable to import data from:%s',filePath);
               end
           end
           
           if ~isempty(sessionSample)
               isSuccessful = 1;
           end
           
           % SEPERATE THE DATA STREAMS
           clearvars -except sessionSample isSuccessful
       end
       % IMPORT A SESSION FROM CYCLE NUMBER
       function [sessionSample] = importCycleData(obj,cycle)
           % This function extracts the DATA, META, EVENTS structures and
           % return them collectively as a sample object.
           
           % INPUT HANDLING
           if ~exist('cycle','var')
               cycle = 1;
           end
           
           % FETCH SESSION FOLDER SET
           illegalPaths = {'.','..'};
           rawDirectories = dir(obj.sessionPath);
           % OMIT KNOWN BAD PATHS
           sessionIndex = rawDirectories(~ismember({rawDirectories.name},illegalPaths));
           % OMIT NON-FOLDERS ( JUST THE LIST OF SESSION DIRECTORIES )
           sessionIndex = sessionIndex([sessionIndex.isdir] == 1);
           
           if isempty(sessionIndex(cycle).name)
               warning('No cycle data available.');
               sessionSample = [];
               return
           end
           
           % CYCLES/SESSIONS ARE ASSUMED IN FOLDER ORDER
           cyclePath = strcat(obj.sessionPath,'\',sessionIndex(cycle).name);
           % IMPORT THE SESSION .MAT FILES
           sessionSample = cell2struct(cell(size(obj.OMASFileNames)),obj.OMASFileNames,2);
           for OMASfile = 1:numel(obj.OMASFileNames)
               filePath = strcat(cyclePath,'\',obj.OMASFileNames{OMASfile},'.mat');
               try
                   load(filePath,'-mat');
                   % IMPORT FILE AS VARIABLE
                   sessionSample.(obj.OMASFileNames{OMASfile}) = eval(obj.OMASFileNames{OMASfile});
               catch
                   warning('Unable to import data from:%s',filePath);
               end
           end
           % SEPERATE THE DATA STREAMS
           clearvars -except sessionSample
       end
       % GET THE COMPLETE SESSION DATA STRUCTURES 
       function [allSamples] = importSession(obj,monteCarloPath)
           % This function examines a generated session directory to
           % extract data for data analysis and plot generation.
           
           % INPUT HANDLING
           if ~exist('monteCarloPath','var')
              monteCarloPath = obj.outputPath; 
           end
           
           % FETCH SESSION FOLDER SET
           illegalPaths = {'.','..'};
           rawDirectories = dir(monteCarloPath);
           totalSessions = numel(rawDirectories) - numel(illegalPaths);
           % OMIT KNOWN BAD PATHS
           sessionIndex = rawDirectories(~ismember({rawDirectories.name},illegalPaths));
           
           % BEGIN MOVING THROUGH THE VALID SESSION DIRECTORIES
           allSamples = struct('DATA',[],'META',[],'EVENTS',[]);
           for sessionID = 1:totalSessions
               % PATH TO EACH OMAS SESSION DIRECTORY
               tempSessionPath = strcat(monteCarloPath,'\',sessionIndex(sessionID).name);
               
               % IMPORT THE SESSION .MAT FILES
               for OMASfile = 1:numel(obj.OMASFileNames)
                   dataPath = strcat(tempSessionPath,'\',obj.OMASFileNames{OMASfile},'.mat');
                   try
                       load(dataPath,'-mat');
                       % IMPORT FILE AS VARIABLE
                       allSamples(sessionID).(obj.OMASFileNames{OMASfile}) = eval(obj.OMASFileNames{OMASfile});
                   catch
                       warning('Unable to import data from:%s',dataPath);
                   end
               end
           end
           
           % SEPERATE THE DATA STREAMS
           clearvars -except allSamples
       end
   end
   % //////////////////// SAMPLE DATA INTERACTIONS ////////////////////////
   methods (Static)
       % GET TRAJECTORY STATISTICS
       function [minSeparation,maxSeparation] = getTrajectoryStatistics(sample)
           % This function is designed to recover relevant trajectory
           % statistics and return them as representative samples.
%            DATA = sample.DATA;
%            SIM  = sample.META;
           
           % BUILD OBJECT SET FOR ALL NON-WAYPOINTS
           objectMETA = sample.META.OBJECTS([sample.META.OBJECTS.type] ~= OMAS_objectType.waypoint);
           [LIA,LOCB] = ismember(sample.META.globalIDvector,[objectMETA.objectID]);
           collidableIDs = LOCB(LIA);
           
           % OUTPUT CONTAINER
           separationTimeSeries = inf(numel(collidableIDs),numel(collidableIDs),sample.META.TIME.numSteps);
           for IDnumA = 1:numel(collidableIDs)
               % THE AGENT STATES
               [objectStatesA] = OMAS_getTrajectoryData(sample.DATA,collidableIDs(IDnumA));
               
               for IDnumB = 1:numel(collidableIDs)
                   if IDnumA == IDnumB
                       % OMIT SEPERATIONS BETWEEN ITSELF
                       continue
                   else
                       % GET THE AGENT STATE TIMESERIES
                       objectStatesB = OMAS_getTrajectoryData(sample.DATA,collidableIDs(IDnumB));
                       centroidSeparations = objectStatesB(1:3,:) - objectStatesA(1:3,:);  % seperation of the centroids
                       centroidSeparations = sqrt(sum(centroidSeparations.^2,1));
                       % STORE IN SEPERATION TIMESERIES MATRIX
                       separationTimeSeries(IDnumB,IDnumA,:) = centroidSeparations;
                   end
               end
           end
           % MIN AND MAXIMUM SEPERATION MATRICES [IDA by IDB]
           %                 minABAxes = collidableIDs;
           minABMatrix = min(separationTimeSeries,[],3);              % Minimum agent-object separations over the timeseries
           maxABMatrix = max(separationTimeSeries,[],3);              % Maximum agent-object separations over the timeseries
           % GET THE MINIMUM SEPARATION FOR ALL COLLIDABLES FOR THE COMPLETE TIMESERIES
           minSeparation = min(minABMatrix(~isinf(minABMatrix)));
           maxSeparation = max(maxABMatrix(~isinf(maxABMatrix)));
       end
       % GET AGENT COMPUTATION STATISTICS 
       function [mean_dt,min_dt,max_dt,dt_timeSeries,timeVector] = getAgentTimingStatistics(sample)
           % This function computes the agent statistics across all
           % Monte-Carlo cycles.
           
           % // CALCULATE THE AGENT COMPUTATION MEANS FOR EACH CYCLE //
           metaObjects = sample.META.OBJECTS;
           agentMetaObjects = metaObjects([metaObjects.type] == OMAS_objectType.agent);
           agentLogicals = ismember(sample.META.globalIDvector,[agentMetaObjects.objectID]);
           % GET THE AGENT-OBJECT SUBSET
           agentIndex  = sample.DATA.objectIndex(agentLogicals);
           agentNumber = sample.META.totalAgents;
           cycleTimeVector = sample.META.TIME.timeVector;
           
           % OUTPUTS
           mean_dt = 0;
           max_dt  = 0;
           min_dt  = 0;
           dt_timeSeries = zeros(size(cycleTimeVector));
           
           % CALCULATE THE AGENT MEANS
           for agentNum = 1:agentNumber
               % DATA LOCAL TO A GIVEN AGENT OF A GIVEN CYCLE
               agentDATA = agentIndex{agentNum}.DATA;
               % THE MEAN NOMINAL COMPUTATION TIME
               mean_dt = mean_dt + agentDATA.algorithm_mean_dt/agentNumber;
               % THE MEAN MAXIMUM COMPUTATION TIME
               max_dt  = max_dt + agentDATA.algorithm_max_dt/agentNumber;
               % THE MEAN MINIMUM COMPUTATION TIME
               min_dt  = min_dt + agentDATA.algorithm_min_dt/agentNumber;
               
               % BUILD CYCLE MEAN COMPUTATION TIME-TIMESERIES
               % The simulations are of varying lengths, therefore the
               % computation times must be padded to the maximum
               % duration to ensure regular data.
               % PAD THE ALGORITHM TIMERS TO THE FULL TIMEVECTOR LENGTH
               %                    algIndicator = zeros(size(cycleTimeVector));
               %                    algIndicator(1:numel(agentDATA.algorithm_indicator)) = agentDATA.algorithm_indicator;
               algCompTime = zeros(size(cycleTimeVector));
               algCompTime(1:numel(agentDATA.algorithm_indicator)) = agentDATA.algorithm_dt;
               
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
           timeVector = sample.DATA.timeVector;
       end
       % GET COLLISION STATISTICS FOR THE GIVEN CYCLES
       function [collisions,waypoints] = getEventStatistics(sample)
           % This function computes the collision statistics across all
           % Monte-Carlo cycles.
           
           % CALCULATE THE MEAN NUMBER OF COLLISIONS
           collisions = sample.DATA.collisions;
           % CALCULATE THE MEAN NUMBER OF ACHIEVED WAYPOINTS
           waypoints = sample.DATA.waypointsAchieved;
       end
   end
   % //////////////////////// CYCLE FUNCTIONS /////////////////////////////
   methods    
       % RUN A DEFINED NUMBER OF PARALLEL CYCLES
       function [obj,cycle,cycleMessage] = parallelCycles(obj)
           % This function runs the specified number of cycles in parallel
           % utilising the Matlab Parallel Toolbox.
           
           fprintf('[%s]\tConfirming parallel worker pool...\n',obj.phase);
           
           % CHECK THE WORKER POOL IS INITIALISED
           obj.threadpoolInterface('RENEW');
           
           fprintf('\n[%s]\tEXECUTING %s PARALLEL MONTE-CARLO CYCLES...\n\n',obj.phase,num2str(obj.cycles));
           cycleMessage = cell(1,obj.cycles);
           parfor cycle = 1:obj.cycles
               fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',obj.phase,cycle);
               try
                   obj.singleCycle();
               catch cycleError
                   warning('\n[%s]\t ...CYCLE %d FAILED.\n\n',obj.phase,cycle);
                   cycleMessage{cycle} = strcat('\n',cycleError.message);
               end
               fprintf('\n[%s]\tCYCLE %d COMPLETE.\n\n',obj.phase,cycle);
           end
           % PARSE MESSAGE LOG
           messageBinary = cellfun('isempty',cycleMessage) ~= 1;           % Flags if message present
           cycle = find(messageBinary,1, 'first');                         % Get the index of the first error
           if isempty(cycle)
               cycle = obj.cycles;
           end
           cycleMessage = cycleMessage(messageBinary);  % Build list of total errors
           cycleMessage = strjoin(cycleMessage);
           
           % KILL THE CURRENT THREAD POOL (FREE MEMORY?)
           obj.threadpoolInterface('KILL');
           
           clearvars -except obj cycle cycleMessage
           
           fprintf('[%s]\t%d PARALLEL CYCLES COMPLETE.\n',obj.phase,obj.cycles);
       end
       % RUN A DEFINED NUMBER OF SEQUENCIAL CYCLES
       function [obj,cycle,cycleMessage] = sequentialCycles(obj)
           % This function runs the specified number of cycles sequentially
           % and returns the collective data for analysis.
           
           fprintf('\n[%s]\tEXECUTING %d SEQUENTIAL MONTE-CARLO CYCLES...\n\n',obj.phase,obj.cycles);
           
           cycleMessage = [];
           for cycle = 1:obj.cycles
               fprintf('\n[%s]\tCYCLE %d STARTING.\n\n',obj.phase,cycle);
               try
                   obj.singleCycle();
               catch cycleError
                   warning('\n[%s]\t ...CYCLE %d FAILED.\n\n',obj.phase,cycle);
                   cycleMessage = cycleError.message;
                   return
               end
               clearvars -except obj cycle cycleMessage
               fprintf('\n[%s]\tCYCLE %d COMPLETE\n\n',obj.phase,cycle);
           end
           
           clearvars -except obj cycle cycleMessage
           
           fprintf('[%s]\t%d SEQUENTIAL CYCLES COMPLETE.\n',obj.phase,obj.cycles);
       end
       % RUN ONE PARAMETERISED SIMULATION CYCLE
       function [MC_DATA,MC_META] = singleCycle(obj)
           % This function computes a single monte-carlo simulation cycle
           % using the monte-carlo defined simulation parameters.
           
           % ASSUMPTIONS:
           % - The parameters of the monte-carlo object are passed to
           %   OpenMAS directly.
           % - The local threadpool option must be disabled to allow the
           %   threads to be used for parallel simulations.
           
           
           % We must perturb the initial conditions of the simulations to
           % allow variation to occur between cycles.
           
           % CREATE A LOCAL COPY OF THE OBJECT INDEX 
           objectIndex = obj.SETTINGS.objects(obj.sessionNumber,:);
           objectIndex = objectIndex(cellfun('isempty',objectIndex) == 0); % Omit empty cells
           
           % PERTURB THE INITIAL XY POSITIONS
           for objectNum = 1:numel(objectIndex)
               try               
                   % GET THE OBJECTS GLOBAL PROPERTIES
                   globalParams = objectIndex{objectNum}.VIRTUAL;
                   % PERTURBATE THE GLOBAL POSITIONS
                   % globalParams.globalPosition = globalParams.globalPosition + obj.SETTINGS.positionSigma*randn(3,1);    % 3D perturbation
                   globalParams.globalPosition = globalParams.globalPosition + [obj.SETTINGS.positionSigma*randn(2,1);0]; % 2D perturbation
                   objectIndex{objectNum}.VIRTUAL = globalParams;
               catch
                   warning('Unable to perturb global position of object %s',num2str(objectNum));
               end
           end
           
           clearvars -except obj objectIndex
           
           % PASS THE INPUTS INTO THE OMAS INITIALISER
           [MC_DATA,MC_META] = OMAS_initialise(...
           'objects',objectIndex,...
           'duration',obj.SETTINGS.duration,... 
           'dt',obj.SETTINGS.dt,...                                        % Timing parameters
           'figures',obj.SETTINGS.figures,...
           'warningDistance',obj.SETTINGS.warningDistance,...
           'conditionTolerance',obj.SETTINGS.conditionTolerance,...
           'visabilityModifier',obj.SETTINGS.visabilityModifier,...
           'outputPath',obj.sessionPath,...
           'threadPool',logical(false),...                                 % Disable the agent-loop thread pool
           'verbosity',obj.SETTINGS.verbosity,...
           'gui',logical(false),...
           'monteCarloMode',logical(true));                                % Enforce monte-carlo file creation
       
            clearvars -except MC_DATA MC_META
       
           % PAUSE FOR STABILITY
           pause(0.1); % 0.1s delay
       end
   end
   % ////////////////////// MONTE-CARLO FIGURES ///////////////////////////
   methods
       % PLOT THE MEAN COMPUTATION TIME TIME-SERIES 
       function [figureHandle] = get_MC_timeSeries(obj,agentData)
           % This function generates the mean computation time - time
           % series. 
           
           % BRING UP THE FIRST SIMPLE
           [cycleSample] = obj.importCycleData(1);
           exampleMETA = cycleSample.META;
           
           % GET THE COMPLETE FIGURE PROPERTY SET
           [figureProperties] = OMAS_figureProperties(exampleMETA.OBJECTS);
           
           % META FIGURE PROPERTIES
           titleString = sprintf('Monte-Carlo (%s cycle) mean computation times',num2str(obj.cycles));
           figurePath = strcat(obj.sessionPath,'\','computationTimeSeries');
           figureHandle = figure('Name','Monte-Carlo Analysis: Agent mean computation times');
           setappdata(figureHandle, 'SubplotDefaultAxesLocation', [0.08, 0.1, 0.90, 0.88]);
           set(figureHandle,'Position',figureProperties.windowSettings);        % [x y width height]
           set(figureHandle,'Color',figureProperties.figureColor);          % Background colour 

           % GENERATE THE FIGURE
           axesHandle = gca;
           displayData = agentData.meanLoopTimeSeries*100;                 % Convert s to ms
           lineHandle = plot(axesHandle,agentData.meanTimeVector,displayData);
           set(lineHandle,'LineWidth',figureProperties.LineWidth);
       
           % FIGURE PROPERTIES
           title(titleString,'fontweight',figureProperties.fontWeight,'fontsize',figureProperties.titleFontSize,'FontSmoothing','on');
           xlabel(axesHandle,'t (s)','fontweight',figureProperties.fontWeight,'fontSize',figureProperties.axisFontSize,'FontSmoothing','on');
           ylabel(axesHandle,'Computation Time (ms)','fontweight',figureProperties.fontWeight,'fontSize',figureProperties.axisFontSize,'FontSmoothing','on');
           set(axesHandle,'FontSize',figureProperties.axisFontSize,'fontWeight',figureProperties.fontWeight,'FontSmoothing','on');
           set(axesHandle,'GridAlpha',0.25,'GridColor','k');
           xlim([ 0 exampleMETA.TIME.endTime]);
           grid on; box on; hold off;
           
           % SAVE FIGURE TO OUTPUT DIRECTORY
           savefig(figureHandle,figurePath);
           
           % SAVE AS PDF
           set(figureHandle,'Units','Inches');
           pos = get(figureHandle,'Position');
           set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
           print(figureHandle,figurePath,'-dpdf','-r0');
       end
       % MONTE-CARLO SIMULATION
       function writeSessionDescriptor(obj) 
           % This function aims to write reference information to the study
           % directory for quick reference.
           
           % WRITE THE OMAS CONFIGURATION TO THE CSV 
           simParams = rmfield(obj.SETTINGS,'objects');
           filePath = strcat(obj.sessionPath,'\OMAS_configuration.csv');
           writetable(struct2table(simParams),filePath);
           
           % WRITE THE OBJECT SET TO THE CSV
           objectIndex = obj.SETTINGS.objects(obj.sessionNumber,:);
           objectIndex = objectIndex(cellfun('isempty',objectIndex) == 0); % Remove empty cells
           
           filePath = strcat(obj.sessionPath,'\OMAS_objectNames.csv');
           fid = fopen(filePath, 'w');
           fprintf(fid,'ObjectID,Name,Type\n');
           for row = 1:numel(objectIndex)
               fprintf(fid,'%d ,%s ,%s\n',(objectIndex{row}.objectID),char(objectIndex{row}.name),char(class(objectIndex{row}))); 
           end
           fclose(fid);
       end
   end
   % //////////////////////// SETUP FUNCTIONS /////////////////////////////
   methods (Static)
       % WRITE LOG-FILE
       function writeLogFile(filePath,entryString)
           % This function writes a time and dated entry to the log file.
           
           filePath = strcat(filePath,'\log.log');
           % GET CURRENT DATE AND TIME
           tString = datestr(datetime('now'),' yyyy-mm-dd @ HH-MM-SS');  
           % WRITE ENTRY TO FILE
           fid = fopen(filePath, 'a+');
           fprintf(fid,'[%s]\t%s\n',tString,entryString); 
           fclose(fid); 
       end
       % CONDITION OUTPUT DIRECTORY
       function [sessionIndex]  = getValidSessionPaths(sessionPath)
           % This function is used to infer an analysis parameters given
           % the files in a provided path. It is assumed that the provided
           % path contains an index of session data sets.

           assert(ischar(sessionPath),'Please provide a valid path string.');
           
           % FETCH SESSION FOLDER SET
           illegalPaths = {'.','..'};
           rawDirectories = dir(sessionPath);
           % OMIT KNOWN BAD PATHS
           sessionIndex = rawDirectories(~ismember({rawDirectories.name},illegalPaths));
           % OMIT NON-FOLDERS ( JUST THE LIST OF SESSION DIRECTORIES )
           sessionIndex = sessionIndex([sessionIndex.isdir] == 1); 
           % GET ONLY CORRECTLY LABELLED SESSION DIRECTORIES
           sessionIndex = sessionIndex(strncmpi({sessionIndex.name},'sessiondata',7));
           % CONVERT TO lIST OF PATHS
           sessionIndex = {sessionIndex.name};
       end
       % BUILD THE RELATIVE MONTE-CARLO DIRECTORY
       function [monteCarloPath] = getMonteCarloDirectory(outputPath,numObjects,numCycles)
           % This function build the monte-carlo data output directory.
           
           % CONFIRM PATH CONVENTION
           if ~strcmp(outputPath(end),'\')
               outputPath = strcat(outputPath,'\');
           end
           
           % GENERATE THE MONTE-CARLO SUBDIRECTORY
           t = datestr(datetime('now'),' yyyy-mm-dd @ HH-MM-SS');          % Record current time
           fileString = sprintf('monteCarlo sessiondata %s %dobjects %dcycles',t,numObjects,numCycles); % Build filestring
           monteCarloPath = strcat(outputPath,fileString);
           % CREATE OUTPUT DIRECTORY
           assert(mkdir(monteCarloPath) == 1,'Unable to create Monte-Carlo output directory');
       end
       
       % CONFIGURE THE PARALLEL/MULTITHREADED WORKERS
       function [poolObject] = threadpoolInterface(conf)
           % This function is designed to generate a threadpool for the forth comming
           % application if one does not already exist. This is only necessary if the
           % aagent loops are sufficiently complicated (i.e. dt > 1s)
           % OUTPUTS:
           % poolObj - The pool object
           
           % PREPARE PARALLEL POOL
           poolObject = gcp('nocreate'); % Get current pool status 
           
           switch upper(char(conf))
               case 'CONFIRM'
                    return
               case 'RENEW'
                   % TERMINATE THE POOL OBJECT
                   delete(poolObject);
                   fprintf('[MONTE-CARLO]\tRenewing parallel pool...\n');
                   % Open new pool for the monte-carlo simulations
                   simCluster = parcluster('local');
                   % CREATE THE PARALLEL POOL OBJECT
                   poolObject = parpool(simCluster,...
                                        'IdleTimeout',30,...               % Generate thread pool, 2hour timeout
                                        'SpmdEnabled',true);
                                    
                   fprintf('[MONTE-CARLO]\t... Done.\n');                 
               case 'KILL'
                   fprintf('[MONTE-CARLO]\tDeleting current active pool...\n');
                   delete(poolObject);
                   fprintf('[MONTE-CARLO]\t... Done.\n');
               otherwise
                   warning('Unrecognised threadpool interface command.');
                   return
           end
           
           fprintf('[MONTE-CARLO]\t... parallel pool ready.\n');
           
%            % PREPARE PARALLEL POOL
%            p = gcp('nocreate'); % Get current pool status
%            if ~isempty(p)
%                % A thread pool exists
%                fprintf('[MONTE-CARLO]\t... Existing parallel pool found.\n');
%                poolObject = p;
%                return
%            else
%                fprintf('[MONTE-CARLO]\t... no active pool found.\n');
%            end
%            
%            % Open new pool for the monte-carlo simulations
%            fprintf('[MONTE-CARLO]\tOpening new parallel pool...\n');
%            simCluster = parcluster('local');
%            poolObject = parpool(simCluster,...
%                                 'IdleTimeout', 30,...   % Generate thread pool, 2hour timeout
%                                 'SpmdEnabled',true);
           
%            fprintf('[MONTE-CARLO]\t... parallel pool ready.\n');
       end
       % PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
       function [config] = configurationParser(defaultConfig,inputParameters)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called in the class constructor
             
            % MOVE THROUGH THE PARAMETER PAIRS (i.e. ,'agents',agentIndex,)
            for parameterIndex = 1:numel(inputParameters)
                % FOR EACH USER-DEFINED PARAMETER
                givenParameter = inputParameters{parameterIndex};
                if ~ischar(givenParameter)
                    continue
                else
                    % IF THE OBJECT HAS A PROPERTY BY THAT NAME
                    if isstruct(defaultConfig) && isfield(defaultConfig,givenParameter)
                        defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};  % Make a substitution
                    end
                    if ~isstruct(defaultConfig) && isprop(defaultConfig,givenParameter)
                        defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};   % Make a substitution
                    end
                end
            end 
            % ERROR CHECKING //////////////////////////////////////////////
            % CHECK AGENTS HAVE BEEN PROVIDED
            if isfield(defaultConfig,'objects') && isempty(defaultConfig.objects)
                error('You must specify and object cell array using the (objects) input parameter'); 
            end
            config = defaultConfig;
       end
       % GET DEFAULT MONTE-CARLO SETTINGS
       function [defaultConfig] = getDefaultConfig()
           % DEFAULT MONTE-CARLO SETTINGS
           defaultConfig = struct(...
               'verbosity',0,...                                           % Append simulation parameters
               'gui',logical(false),...
               'figures',{{'none'}},...
               'objects',[],...
               'warningDistance',10,...
               'conditionTolerance',1E-3,...
               'visabilityModifier',1,...
               'positionSigma',0,...
               'startTime',0,...
               'endTime',inf,...
               'duration',10,...
               'dt',0.1);
           % RETURN THE INITIALISED SETTINGS
       end
   end
end