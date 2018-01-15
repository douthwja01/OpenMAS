%% OPENMAS SCENARIO SIMULATOR (OMAS_initialise.m) %%%%%%%%%%%%%%%%%%%%%%%%%
% This function is designed to accept a set of agent objects to inject into
% the simulation instance. The program will then iterate through the
% scenario time steps.

% Author: James A. Douthwaite 16/05/2016

% MAIN SIMULATION SCHEDULE
function [ DATA,META ] = OMAS_initialise(varargin)
% This function handles conducts the precalculations and setup functions
% necessary for the simulation to be ran on the agentIndex
% INPUTS:
% SIM         - The simulation configuration structure
% objectIndex - The vector of object classes

% OUTPUTS:
% DATA       - The simulation output data structure

% PARSE THE INPUT STRUCTURE
[SIM,objectIndex] = validateSimulationInputs(varargin);
if ~isstruct(SIM)
    warning('Simulation input parameters are invalid.');
    DATA = 0;
    return
end

%% SIMULATION CONFIGURATION PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[%s]\tINITIALISING SIMULATION...\n\n',SIM.phase);
close all; % <<< close all figures

% GET THE SIMULATIONS SUBDIRECTORIES (if not already added to path)
directorySet = {'environment','toolboxes','events','objects','scenarios'};
primaryDir = pwd; 
getSimulationDirectories(primaryDir,directorySet);

% MAKE ANY ADDITIONAL PREPERATIONS TO THE META DATA
clear META
[META,objectIndex] = defineSimulationParameters(SIM,objectIndex);
clearvars SIM;

% CONFIGURE THE OUTPUT SETTINGS
fprintf('[%s]\tOUTPUT CONFIGURATION:\n',META.phase);
% DEFINE THE RELATIVE OUTPUT PATH
relativeOutputDir = ['\data\',META.outputPath]; % Relative data directory
try
    [META.outputPath,META.outputFile] = configureOutputDirectory(primaryDir,relativeOutputDir);  % Get current working working directory variables
    fprintf('[%s]\tSession output directory created: %s.\n',META.phase,META.outputFile);
    fprintf('[%s]\tOutput path: %s\n',META.phase,META.outputPath); 
catch
    warning('[%s}\tFailed to generate output directory: %s',META.phase,META.outputPath);
    DATA = 0;
    SIM = META;
    return
end
% SAVE THE META DATA OBJECT
save([META.outputPath,'META.mat'],'META'); % Record the simulation parameters

% INITIALISE THE VR-WORLD IF REQUESTED
if META.visualisation
    initialiseVisualisation()
end

% DISPLAY SIMULATION CONFIGURATION 
fprintf('[%s]\tCONFIGURATION SUMMARY:\n',META.phase);
fprintf('[%s]\tObjects: %s\tAgents: %s\tObstacles: %s\tWaypoints: %s\n',...
        META.phase,num2str(META.totalObjects),num2str(META.totalAgents),...
        num2str(META.totalObstacles),num2str(META.totalWaypoints));                    % Object summary
fprintf('[%s]\tDuration: %ss\tSampling: %ss\tSteps: %s\n',META.phase,...
        num2str(META.TIME.simTime),num2str(META.TIME.dt),num2str(META.TIME.numSteps)); % Timing summary
   
%% //////// BEGIN RUNNING THROUGH THE SIMULATION TIME CYCLE ///////////////
fprintf('[%s]\tMoving to simulation...\n',META.phase);
[DATA,META,EVENTS,objectIndex] = OMAS_mainCycle(META,objectIndex);

% BEGIN RUNNING EXTERNAL DATA ANALYSIS
try
    fprintf('[%s]\tMoving to post-simulation analysis...\n[%s]\n',META.phase,META.phase);
    [DATA] = OMAS_analysis(META,objectIndex,EVENTS,DATA);                  % Jump to external analysis program
catch analysisError
    warning('[ERROR] A problem occurred interpreting the output data');
    warning(analysisError.message);
end
fprintf('\n[%s]\tExiting...\n\n',META.phase);
end

%% SIMULATION SETUP FUNCTIONS /////////////////////////////////////////////
% VALIDATE THE INPUT STRUCTURE & DEFAULT PARAMETERS
function [ SIM,objectIndex ] = validateSimulationInputs(inputSet)
% This function is to accept the user provided input variable set and parse
% any specified parameter and apply it to the SIM structure. 
% INPUTS:
% inputSet    - The cell array of simulation inputs, in string label, value
%               pairs.
% OUTPUTS:
% SIM         - The complete input, with default values
% objectIndex - The cell array of simulation objects

%% DEFINE DEFAULT INPUT CONDITIONS
% DEFAULT SIMULATION CONDITIONS
objectIndex = {};
defaultFigureSet = {'EVENTOVERVIEW','FIG'};
SIM = struct('warningDistance',10,...                                      % Default to 10m seperation
            'visabilityFactor',100,...                                     % Visability percentage
                     'figures',{defaultFigureSet},...                      % Output default figure set
                  'outputPath','',...                                      % Default to the top level of /data
                  'threadPool',boolean(0),...                              % Disable threadpool by default
                 'verboseMode',boolean(1),...                              % Default to displaying debug information
               'visualisation',boolean(0));                                % Output visualisation
             
% DEFAULT TIMING PARAMETERS
SIM.TIME = struct('simTime',10,'dt',0.1,'startTime',0);                    % Assign default timing parameters        

%% BEGIN PARSING INPUT VARIABLES  /////////////////////////////////////////         
fprintf('[SETUP]\tConfirming input variables.\n');

%% TIMING PARAMETERS //////////////////////////////////////////////////////
tmp = strncmpi(inputSet,'SIMTIME',7)|strncmpi(inputSet,'RUNTIME',7); 
if any(tmp)
    SIM.TIME.simTime = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.TIME.simTime) == 1,'The simulation total time "simTime" must be numeric.');
end
tmp = strncmpi(inputSet,'DT',2)|strcmp(inputSet,'TIMESTEP'); 
if any(tmp)
    SIM.TIME.dt = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.TIME.dt) == 1,'The simulation time step "dt" must be numeric.');
end
tmp = strcmpi(inputSet,'STARTTIME')|strncmpi(inputSet,'T0',2); 
if any(tmp)
    SIM.TIME.simTime = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.TIME.simTime) == 1,'The simulation total time "startTime" must be numeric.');
end

%% THE OBJECT/AGENT SET ///////////////////////////////////////////////////
tmp = strcmpi(inputSet,'OBJECTS')|strncmpi(inputSet,'AGENTS',6);  
if any(tmp)
    tmp = inputSet{find(tmp) + 1};
    assert(iscell(tmp) == 1,'The object set must be specified as a cell column vector.');
    objectIndex = vertcat(objectIndex,tmp);
end

%% OUTPUT DATA ////////////////////////////////////////////////////////////
% PLOT REQUESTS
tmp = strncmpi(inputSet,'FIGURES',7)|strncmpi(inputSet,'PLOTS',5); 
if any(tmp)
    SIM.figures = inputSet{find(tmp) + 1}; 
    assert(iscell(SIM.figures) == 1,'Figure set must be specified as a cell array of string names.');
end
% OUTPUT DIRECTORY/PATH
tmp = strncmpi(inputSet,'DIRECTORY',5)|strncmpi(inputSet,'OUTPUTPATH',6);  
if any(tmp)
    SIM.outputPath = inputSet{find(tmp) + 1};
    assert(ischar(SIM.outputPath) == 1,'Please specify a relative directory string.');
end

%% GENERAL CONFIGURATION //////////////////////////////////////////////////
% GLOBAL OBJECT SEPERATION WARNING RADIUS
tmp = strncmpi(inputSet,'WARNING',4)|strncmpi(inputSet,'NEARMISS',7); 
if any(tmp)
    SIM.warningDistance = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.warningDistance) == 1,'Specify the warning distance as a numeric (in meters).');
end
% VISABILTY REDUCTION RATIO (IN FUTURE PROPORTIONAL TO ISO PROPERTIES)
tmp = strcmpi(inputSet,'VISABILITY'); 
if any(tmp)
    SIM.visabilityFactor = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.visabilityFactor) == 1,'Specify the visability factor as a numeric ratio.');
end
% DEFINE THE SIMUATION OUTPUT MODE
tmp = strncmpi(inputSet,'VERBOSE',6)|strncmpi(inputSet,'SILENT',3); 
if any(tmp)
    tmp = inputSet{find(tmp) + 1};
    assert(isnumeric(tmp) == 1,'Verbosity must be specified as a logical.');
    if tmp > 0
        SIM.verboseMode = boolean(1); % Assume a value higher than one is a request
    else
        SIM.verboseMode = boolean(0);
    end   
end
% MULTI-THREADING TASK HANDLING
tmp = strcmpi(inputSet,'THREADPOOL')|strncmpi(inputSet,'THREADING',4);
if any(tmp)
    tmp = inputSet{find(tmp) + 1};
    assert(isnumeric(tmp) == 1,'Multi-threading is boolean enabled (0:1).');
    if tmp > 0
        SIM.threadPool = boolean(1); % Assume a value higher than one is a request
    else
        SIM.threadPool = boolean(0);
    end     
end

% VISUALISATION MODE
tmp = strcmpi(inputSet,'VISUALISATION')|strncmpi(inputSet,'VISUALISE',5);
if any(tmp)
    tmp = inputSet{find(tmp) + 1};
    assert(isnumeric(tmp) == 1,'Visualisation mode is boolean enabled (0:1).');
    if tmp > 0
        SIM.visualisation = boolean(1); % Assume a value higher than one is a request
    else
        SIM.visualisation = boolean(0);
    end     
end

% SIMULATION PHASE INDICATOR
SIM.phase = 'SETUP'; % Initialise in setup mode
fprintf('[%s] --> Input parameters defined.\n',SIM.phase);
end
% PREPARE THE META STRUCTURE AND INITIALISE OBJECT INDEX
function [ SIM,objectIndex ] = defineSimulationParameters(SIM,objectIndex)
% INPUTS:
% SIM         - The simulation inputs
% objectIndex - A vector of object classes
% OUTPUTS:
% SIM - The prepared simulations parameters

%% GENERAL PARAMETER INITIALISATION
% DETERMINE IF A PARALLEL POOL IS NEEDED
if SIM.threadPool 
    [SIM.threadPool] = configureThreadpool();
end

%% SIMULATION PARAMETER INITIALISATION
% DEFINE TIME PARAMETERS
SIM.TIME.endTime    = SIM.TIME.startTime + SIM.TIME.simTime;               % Define the end time
SIM.TIME.timeVector = SIM.TIME.startTime:SIM.TIME.dt:(SIM.TIME.endTime);   % The progressive time vector
SIM.TIME.numSteps   = size(SIM.TIME.timeVector,2);                         % Number of resulting steps
SIM.TIME.frequency  = 1/SIM.TIME.dt;                                       % Define the simulation frequency 

% INSTANTEOUS PARAMETERS
SIM.TIME.currentTime = 0;
SIM.TIME.currentStep = 1;

% DEFINE THE AVAILABLE STATUS TYPES FROM THE 'eventType' ENUMERATIONS
[K,SIM.eventNames] = enumeration('eventType');

SIM.eventNames = SIM.eventNames(2:end);
SIM.eventTypes = length(SIM.eventNames);                                   % Size the number of available sim events (neglect base definition)

% DEFINE THE NUMBER OF AGENTS
SIM.totalObjects = numel(objectIndex);                                      % Define the agent number
SIM.totalAgents = 0;
SIM.totalObstacles = 0;
SIM.totalWaypoints = 0;

% DEFINE LOCAL STATE INDICES
positionIndices = 1:3;
rotationIndices = 4:6;
velocityIndices = 7:9;
omegaIndices = 10:12;

%PREPARE THE META.OBJECT GLOBAL STATUS VECTOR
for i = 1:SIM.totalObjects
    % THE OBJECTS POSSESS INITIAL GLOBAL PROPERTIES (defined by
    % the scenario object) THESE ARE:
    % .position     - Global position in East-North-Up (ENU) coordinates
    % .velocity     - Global velocity in ENU coordinates
    % .quaternion   - Global (unit triad) to rotated body axis frame
    
    %% EXTRACT THE GLOBAL POSITION AND VELOCITY
    % The scenario builder provides the position, velocity and quaternion
    % in the ENU frame.
    globalXYZPosition = objectIndex{i}.globalPosition;                     % Get ENU global position
    globalXYZVelocity = objectIndex{i}.globalVelocity;                     % Get ENU global velocity
    staticBodyToRotatedGlobalQuaternion = objectIndex{i}.quaternion;       % REDEFINE INITIAL QUATERNION AS STATIC BODY TO ROTATED EARTH QUATERNION
    
    % CONFIRM INPUTS
    assert(numel(globalXYZPosition) == 3,'Incorrect global position dimensions.');
    assert(numel(globalXYZVelocity) == 3,'Incorrect global velocity dimensions.');
    assert(numel(staticBodyToRotatedGlobalQuaternion) == 4,'Incorrect quaternion dimensions.');
    
    % GET THE INITIAL ROTATIONS
    [R_BG,R_GB] = OMAS_axisTools.quaternionToRotationMatrix(staticBodyToRotatedGlobalQuaternion);
    % CALCULATE THE EULER HEADING
    XYZrotations = rotm2eul(R_BG,'ZYZ');
    
    % ASSIGN INITIAL LOCAL FORWARD-RIGHT-DOWN STATE VECTOR
    localFLUState = zeros(12,1);
    localFLUState(rotationIndices) = XYZrotations;                         % Get the initial rotations with respect to the global axis
    localFLUState(velocityIndices) = R_GB*globalXYZVelocity;               % Get the initial local velocity vector
          
    % CONVERT THE FLU (LOCAL ENU) STATES TO THE FRD (LOCAL NED) 
    localFRDState = localFLUState;
    localFRDState(positionIndices(2)) = -localFRDState(positionIndices(2)); % + left in FLU -> -ve in FRD  
    localFRDState(positionIndices(3)) = -localFRDState(positionIndices(3)); % + up in FLU -> -ve in FRD 
    localFRDState(velocityIndices(2)) = -localFRDState(velocityIndices(2)); % ""
    localFRDState(velocityIndices(3)) = -localFRDState(velocityIndices(3)); % ""
    
    % CONVERT THE LOCAL ANGULAR (FLU) STATES TO FRD (LOCAL NED)
    localFRDState(omegaIndices(2)) = -localFLUState(omegaIndices(2));      % Pitch rotation correction
    localFRDState(omegaIndices(3)) = -localFLUState(omegaIndices(3));      % Yaw rotation correction
    % ASSIGN INITIAL LOCAL NORTH-EAST-DOWN STATE
    objectIndex{i}.localState = localFRDState; 
    
    % IF THERE IS NO DETECTION RANGE GIVEN (object/blind agent)
    if ~isfield(objectIndex{i}.VIRTUAL,'detectionRange')
        objectIndex{i}.VIRTUAL.detectionRange = 0;
    end
    % OBJECTS is a vector of structures recording the cuurent status of each object 
    SIM.OBJECTS(i) = struct('objectID',objectIndex{i}.objectID,...                  % ID reference number
                                'name',objectIndex{i}.name,...                      % Name reference
                               'class',class(objectIndex{i}),...                    % Class reference
                                'type',objectIndex{i}.VIRTUAL.type,...              % Object type
                              'colour',objectIndex{i}.VIRTUAL.colour,...            % Object virtual colour
                              'symbol',objectIndex{i}.VIRTUAL.symbol,...            % Representative symbol (from type)
                                'size',objectIndex{i}.VIRTUAL.size,...              % The objects critical radius
                      'detectionRange',objectIndex{i}.VIRTUAL.detectionRange,...    % The simulations detection horizon for the aircraft
                         'globalState',[globalXYZPosition;globalXYZVelocity;staticBodyToRotatedGlobalQuaternion],...
                                'R_GB',R_GB,...                                     % Current Global-body rotation matrix
                                'R_BG',R_BG,...                                     % Current Body-global rotation matrix
                   'relativePositions',zeros((SIM.totalObjects),3),...              % Relative distances between objects [dx;dy;dz]*objectCount
                 'euclideanSeperation',zeros((SIM.totalObjects),1),...              % Euclidean scalar distances between objects (excluding itself)
                        'objectStatus',zeros((SIM.totalObjects),(SIM.eventTypes/2))); % A matrix recording all interactions against the object

    % ASSIGN NaN TO LOCATIONS OF SELF-REFERENCE            
    SIM.OBJECTS(i).relativePositions((objectIndex{i}.objectID),:) = NaN;
    SIM.OBJECTS(i).euclideanSeperation((objectIndex{i}.objectID),1) = NaN;
    SIM.OBJECTS(i).objectStatus((objectIndex{i}.objectID),:) = NaN;
    
    % COUNT NUMBER OF OBJECT TYPES 
    if SIM.OBJECTS(i).type == OMAS_objectType.agent
        SIM.totalAgents = SIM.totalAgents + 1;
    elseif SIM.OBJECTS(i).type == OMAS_objectType.obstacle
        SIM.totalObstacles = SIM.totalObstacles + 1;
    elseif SIM.OBJECTS(i).type == OMAS_objectType.waypoint
        SIM.totalWaypoints = SIM.totalWaypoints + 1;
    else
        SIM.totalOther = SIM.totalWaypoints + 1; 
    end
    
    % CHECK CONTROL FREQUENCY IS VALID
    % The user specifies an agent control frequency, which must be lower
    % than the minimum simulation sample frequency (1/dt)
    if isprop(objectIndex{i},'sampleFrequency')
        if isinf(objectIndex{i}.sampleFrequency)
            objectIndex{i}.sampleFrequency = (1/SIM.TIME.dt);
        elseif objectIndex{i}.sampleFrequency > (1/SIM.TIME.dt)
            error('Agent %s control frequency (%sHz) cannot be greater than the simulation frequency (%sHz).',...
                  objectIndex{i}.name,num2str(objectIndex{i}.sampleFrequency),num2str(1/SIM.TIME.dt));
        end
    end
end 

% GET THE OBJECT ID REFERENCE VECTOR
SIM.globalIDvector = [SIM.OBJECTS.objectID];
end

%% PRE-SIMULATION TOOLS ///////////////////////////////////////////////////
% INITIALISE THE VR TEST WORLD IF REQUIRED
function initialiseVisualisation(fileName)
    % This function is designed to initialise the generic VR test world
    % file to allow the scene to be updated during the simulation runtime.
    
    
end

% GET THE SIMULATION SUBDIRECTORY SET
function getSimulationDirectories(wrkDir,directorySet)
    % Check the directory set is valid
    if ~iscell(directorySet)
        warning('Simulation subpath vector is invalid.');
        return
    end   
    
    % Move through the sub-directory list
    for i = 1:size(directorySet,2)
        pathName = strcat(wrkDir,'\',directorySet{i});
        try
            addpath(pathName);
        catch simDirectoryError
            warning('Unable to find critical sub-directories, have you changed directories?');
            error(simDirectoryError.message);
        end
    end
end
% CONFIGURE THE PARALLEL/MULTITHREADED WORKERS
function [poolObject] = configureThreadpool()
% This function is designed to generate a threadpool for the forth comming
% application if one does not already exist. This is only necessary if the
% aagent loops are sufficiently complicated (i.e. dt > 1s)
% OUTPUTS:
% poolObj - The pool object

% PREPARE PARALLEL POOL
fprintf('[SETUP]\tConfiguring thread pool...\n');
p = gcp('nocreate'); % Get current pool status
if ~isempty(p)  
   % A thread pool exists
   fprintf('[SETUP]\t--> Existing thread pool found.\n'); 
   poolObject = p;
   return
else
   fprintf('[SETUP]\t--> no active pool found.\n');
end

% Open new pool for the monte-carlo simulations
fprintf('[SETUP]\tOpening new thread pool...\n');
simCluster = parcluster('local');
poolObject = parpool(simCluster,'IdleTimeout', 120);   % Generate thread pool, 2hour timeout
fprintf('[SETUP]\t... Thread pool ready.\n');
end
% CONFIGURE THE OUTPUT DIRECTORY
function [outpath,filestr] = configureOutputDirectory(workdir,relOutputDir)
% INPUTS:
% workdir      - The current working directory string
% relOutputDir - The relative output directory 
% OUTPUTS:
% outpath - The current sessions output path
% filestr - The name of the session output file

% Configure the output directory
%savedir = '\monte_carlo\';
% savedir = strcat(workdir,relOutputDir,'\');
savedir = strcat(workdir,relOutputDir);

% Generate a new sub-directory
t = datestr(datetime('now'),' yyyy-mm-dd @ HH-MM-SS');  % Record current time
filestr = strcat('sessiondata',t);                    % Build filestring
outpath = strcat(savedir,filestr,'\');                  % Define the output directory
flag = mkdir(savedir,filestr);
end