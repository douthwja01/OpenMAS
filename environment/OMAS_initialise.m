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

% ////////// GET OMAS DEPENDENCIES (if not already added to path) /////////
[flag,repoPath] = OMAS_system.GetFileDependancies();
if ~flag
    error('Error getting the OpenMAS file dependancies.');
end

% Show logo
type(strcat(repoPath,"environment\assets\logo.txt"));

% ///////////////// ASSIGN DEFAULT SIMULATION PARAMETERS //////////////////
% The user must be able to run a simulation for an intended amount of time,
% or run with a maximal time for all goals to be complete. 

% /////////////////////// DEFAULT TIMING STRUCTURE ////////////////////////
TIME = struct('duration',1E3,...                                           % Assign default simulation duration (for output matrix scaling)
                    'dt',0.5,...
             'startTime',0,...
           'idleTimeOut',1);
[TIME] = GetParameterOverrides(TIME,varargin);  

% //////////////// DEFAULT GENERAL CONFIGURATION STRUCTURE ////////////////
defaultConfig = struct('figures',{{'events','fig'}},...
               'warningDistance',5,...
            'conditionTolerance',1E-3,...
            'visabilityModifier',2,...
                    'outputPath',OMAS_system.GetOSPathString([pwd,'\data']),... 
                    'systemFile','temp.mat',...
                         'phase','SETUP',...                               % Initialise in setup mode
                    'threadPool',false,...
                'monteCarloMode',false,...
                   'publishMode',false,...
                     'verbosity',2,...
                           'gui',true,...
                       'objects',[],...
                          'TIME',TIME);
                 
% ///////////////////// PARSE USER INPUT PARAMETERS ///////////////////////
fprintf('[SETUP]\tConfirming input variables.\n');
[SIM] = GetParameterOverrides_recursive(defaultConfig,varargin);
clear TIME defaultConfig

% Input sanity checks
assert(islogical(SIM.monteCarloMode),'Please specify the monte-carlo mode as a logical.');
assert(islogical(SIM.gui),'Please specify the use of gui elements as a logical.');
assert(ischar(SIM.systemFile),'Please provide a valid file path.');
assert(ischar(SIM.outputPath),'Please provide a valid output path.');

SIM.threadPool = logical(SIM.threadPool);
assert(islogical(SIM.threadPool),'Please specify the multi-threading mode as a logical.');
SIM.verbosity = uint8(SIM.verbosity);
if ~iscell(SIM.figures) 
    SIM.figures = {SIM.figures};
end
 
% /////////////// SOME PARAMETER DEDUCTION/INPUT HANDLING /////////////////
assert(~isempty(SIM.objects),'Please provide a agent/object object vector.');
% TRIM ANY EMPTY OBJECT CONTAINERS
for ind = 1:numel(SIM.objects)
   if ~isobject(SIM.objects{ind})
      SIM.objects(ind) = {''};                                             % Remove invalid classes from the container
   end
end
% REMOVE EMPTY ENTRIES FROM THE OBJECT ARRAY
SIM.objects = SIM.objects(~cellfun('isempty', SIM.objects)); 

% DISABLE AGENT-WORKER ALLOCATION IN MONTE-CARLO MODE
if SIM.monteCarloMode 
    SIM.threadPool = logical(false);
end

% SIMULATION PHASE INDICATOR
fprintf('[%s] --> Input parameters defined.\n',SIM.phase);

% ///////////////// SIMULATION CONFIGURATION PROCEDURE ////////////////////
fprintf('\n[%s]\tINITIALISING SIMULATION...\n\n',SIM.phase);
close all; 

fprintf('[%s]\tOMAS CONFIGURATION:\n',SIM.phase);
clear META
[META,objectIndex] = ConfigureSimulationParameters(SIM);
clearvars SIM;
fprintf('[%s]\tObjects parameterized.\n',META.phase);

% CONFIGURE THE OUTPUT SETTINGS ///////////////////////////////////////////
fprintf('[%s]\tOUTPUT CONFIGURATION:\n',META.phase);
try % Define output path
    [META.outputPath,META.outputFile] = ConfigureOutputDirectory(META.monteCarloMode,META.outputPath);  % Get current working working directory variables
catch outputPathError
    warning('[%s]\tUnable to generate output directory; %s',META.phase,META.outputPath);
    rethrow(outputPathError);
end    
fprintf('[%s]\tSession output directory created: \n[%s]\t"%s" \n',META.phase,META.phase,META.outputFile);
fprintf('[%s]\tOutput path: \n[%s]\t"%s"\n',META.phase,META.phase,META.outputPath); 
% Save the initial META object 
save([META.outputPath,'META.mat'],'META'); % Record the simulation parameters
% PUSH FIGURE REQUEST TO FILE
figureList = META.figures;
save([META.outputPath,META.systemFile],'figureList');
META = rmfield(META,'figures');            % Remove the figure list from the META structure
clear figureList

% DISPLAY SIMULATION CONFIGURATION 
fprintf('[%s]\tCONFIGURATION SUMMARY:\n',META.phase);
fprintf('[%s]\tObjects: %d\tAgents: %d\tObstacles: %d\tWaypoints: %d\n',...
        META.phase,META.totalObjects,META.totalAgents,META.totalObstacles,META.totalWaypoints);                     % Object summary
fprintf('[%s]\tMax Steps: %d\t Max Duration: %.2fs\t dt: %.3fs\t\n',META.phase,...
        META.TIME.numSteps,META.TIME.duration,META.TIME.dt);   % Timing summary
   
%% //////// BEGIN RUNNING THROUGH THE SIMULATION TIME CYCLE ///////////////
fprintf('[%s]\tMoving to simulation...\n',META.phase);
OMAS_start = tic();
[DATA,META,EVENTS,objectIndex] = OMAS_process(META,objectIndex);
OMAS_finish = toc(OMAS_start); 
fprintf('[%s]\tOperation lasted %.2fs.\n',META.phase,OMAS_finish);

%% ///////////// ONCE COMPLETE - RUN EXTERNAL ANALYSIS ////////////////////
try
    fprintf('[%s]\tMoving to post-simulation analysis...\n[%s]\n',META.phase,META.phase);
    [DATA] = OMAS_analysis(META,objectIndex,EVENTS,DATA);                  % Jump to external analysis program
catch analysisError
    warning('[ERROR] A problem occurred interpreting the output data');
    rethrow(analysisError);
end
% DELETE THE SYSTEM TEMPORARY FILE
delete([META.outputPath,META.systemFile]);

fprintf('\n[%s]\tExiting...\n\n',META.phase);
end

%% SIMULATION SETUP FUNCTIONS /////////////////////////////////////////////
% PREPARE THE META STRUCTURE AND INITIALISE OBJECT INDEX
function [SIM,objectIndex] = ConfigureSimulationParameters(SIM)
% INPUTS:
% SIM         - The simulation inputs
% objectIndex - A vector of object classes
% OUTPUTS:
% SIM         - The prepared simulations parameters

% GENERAL PARAMETER INITIALISATION
% DETERMINE IF A PARALLEL POOL IS NEEDED
if SIM.threadPool 
    [SIM.threadPool] = ConfigureThreadpool();
end

% ////////////// SIMULATION PARAMETER INITIALISATION //////////////////////
% DEFINE TIME PARAMETERS
SIM.TIME.endTime     = SIM.TIME.duration + SIM.TIME.startTime;              
SIM.TIME.timeVector  = SIM.TIME.startTime:SIM.TIME.dt:(SIM.TIME.endTime);  % The progressive time vector
SIM.TIME.currentTime = SIM.TIME.startTime;
SIM.TIME.frequency   = 1/SIM.TIME.dt;                                      % Define the simulation frequency 
SIM.TIME.numSteps    = size(SIM.TIME.timeVector,2);                        % Number of resulting steps
SIM.TIME.currentStep = 1;                                                  % Initialise the current time properties 

% INPUT LOGIC VERIFICATION
assert(SIM.TIME.duration  >= SIM.TIME.dt,'Simulation timestep must be lower than the total duration.');
assert(SIM.TIME.startTime <= SIM.TIME.endTime,'Simulation end-time cannot be before its start-time.');
assert(SIM.TIME.idleTimeOut >= 0,'Simulation time-out duration must be a positive scalar.');

% DEFINE THE AVAILABLE STATUS TYPES FROM THE 'eventType' ENUMERATIONS
[eventEnums,~] = enumeration('eventType');
% SEPERATE THE OBJECT VECTOR
objectIndex = SIM.objects; 
SIM = rmfield(SIM,'objects');                                              % Remove cell array to reduce data overhead
% DEFINE THE NUMBER OF AGENTS
SIM.totalObjects = uint16(numel(objectIndex));                             % Define the object number

%PREPARE THE META.OBJECT GLOBAL STATUS VECTOR
for entity = 1:SIM.totalObjects
    % THE OBJECTS ARE ASSIGNED THEIR GLOBAL PARAMETERS IN THE GLOBAL XYZ 
    % (MATLAB) AXIS SYSTEM;
    % .position     - Global position in XYZ coordinates 
    % .velocity     - Global velocity in XYZ coordinates
    % .quaternion   - Quaternion describing the XYZ rotation
    
    % Extract the objects GLOBAL structure for tear-down
    objectGLOBAL = objectIndex{entity}.GetGLOBAL();
    % Confirm the objects global representation
    assert(IsColumn(objectGLOBAL.position,3)  ,'Global position must be given as a vector [3x1].');
    assert(IsColumn(objectGLOBAL.velocity,3)  ,'Global velocity must be given as a vector [3x1].');
    assert(IsColumn(objectGLOBAL.quaternion,4),'Global attitude must be give as a quaternion vector [4x1].');      
    
    %% ////////////////////// ROTATION CONVENTION /////////////////////////
    % INPUTS:
    % globalXYZPosition     - The global cartesian position
    % globalXYZVelocity     - The global cartesian velocity
    % globalXYZquaternion   - The quaternion defining the rotation of the
    %                         body away from alignment with the global
    %                         axes.                      
    
    % USING MATLABS TOOLBOX COP-OUT
    Rxyz              = quat2rotm(objectGLOBAL.quaternion');                   % from fixed to rotated world
%     localXYZrotations = quat2eul(objectGLOBAL.quaternion','XYZ');
    localXYZrotations = quat2eul(objectGLOBAL.quaternion','ZYX');
    localXYZrotations = localXYZrotations';
    localXYZVelocity  = Rxyz*objectGLOBAL.velocity;                            % Move from rotated vector to axis aligned vector
    objectGLOBAL.R    = Rxyz;
    % Agents idle status is reset before running
    if objectGLOBAL.type == OMAS_objectType.agent
       objectGLOBAL.idleStatus = false;                                    % Initialise the agent as active
    end
    % IF THERE IS NO DETECTION RANGE GIVEN (object/blind agent)
    if ~isfield(objectGLOBAL,'detectionRadius')
        objectGLOBAL.detectionRadius = 0;
    end 
    % Re-assign the modified GLOBAL structure
    objectIndex{entity}.SetGLOBAL(objectGLOBAL);
    
    %% Initialise the object (i.e. object.setup())
    try
        objectIndex{entity} = objectIndex{entity}.setup(localXYZVelocity,localXYZrotations);
    catch initialisationError
        warning('Unable to initialise the local state of object %s :%s',objectIndex{entity}.name);
        rethrow(initialisationError)
    end
    
    % CHECK FOR FORBIDDEN CHARACTERS IN OBJECT NAME ASSIGNMENTS 
    % This prevents issues from occurring with filenames and figures.
    assertionString = sprintf('Forbidden character (:*?"<>|) in object identifier: "%s"',objectIndex{entity}.name);
    assert(~any(contains(objectIndex{entity}.name,{'\','/',':','*','?','"','<','>','|'},'IgnoreCase',true)),assertionString);
        
    %% Assemble the global META object container
    % OBJECTS is a vector of structures containing the object META data 
    SIM.OBJECTS(entity) = struct('objectID',objectIndex{entity}.objectID,...    % ID reference number
                                     'name',objectIndex{entity}.name,...        % Name reference
                                    'class',class(objectIndex{entity}),...      % Class reference
                                     'type',uint8(objectGLOBAL.type),...        % Object type
                                   'hitBox',uint8(objectGLOBAL.hitBoxType),...  % Record the hitbox type
                                   'colour',objectGLOBAL.colour,...             % Object virtual colour
                                   'symbol',objectGLOBAL.symbol,...             % Representative symbol (from type)
                                   'radius',objectGLOBAL.radius,...             % The objects critical radius
                          'detectionRadius',objectGLOBAL.detectionRadius,...    % The simulations detection horizon for the aircraft
                               'idleStatus',objectGLOBAL.idleStatus,...         % Agent completed Task Flag
                              'globalState',[objectGLOBAL.position;objectGLOBAL.velocity;objectGLOBAL.quaternion],...
                                        'R',objectGLOBAL.R,...                  % Current Global-Body rotation matrix
                        'relativePositions',zeros(SIM.totalObjects,3),...       % Relative distances between objects [dx;dy;dz]*objectCount
                             'objectStatus',false(SIM.totalObjects,(numel(eventEnums)-1)/2));                   
    % ASSIGN NaN TO LOCATIONS OF SELF-REFERENCE            
    SIM.OBJECTS(entity).relativePositions(entity,:) = NaN;                 % Mark self reference in the META structure as a NaN.
    % CHECK FOR SIMULATION/SAMPLING FREQUENCY ERROR
    % The user specifies an agent control frequency, which must be lower
    % than the minimum simulation sample frequency (1/dt)
    if isprop(objectIndex{entity},'SENSORS') && isfield(objectIndex{entity}.SENSORS,'sampleFrequency')
%         assert(objectIndex{entity}.SENSORS.sampleFrequency  > (1/SIM.TIME.dt),...
%                'Agent %s sample frequency (%sHz) cannot be greater than the simulation frequency (%sHz).');
        % SPECIAL CASE - UNLIMITED SENSING RATE
        if isinf(objectIndex{entity}.SENSORS.sampleFrequency)
        	objectIndex{entity}.SENSORS.sampleFrequency = (1/SIM.TIME.dt);
        end
    end
end

% GET THE TOTALS FOR EACH OBJECT TYPE FOR REFERENCE
[objEnums,objNames] = enumeration('OMAS_objectType');                      % Define the number of other object types
typeIndex = [SIM.OBJECTS.type];
for entity = 1:numel(objEnums)
    typeLabel = strcat('total',char(objNames(entity)),'s');
    typeLabel(6) = upper(typeLabel(6));
    SIM.(typeLabel) = uint8(sum(typeIndex == objEnums(entity)));
end
% GET THE OBJECT ID REFERENCE VECTOR
SIM.globalIDvector = [SIM.OBJECTS.objectID];
% ORGANISE THE FIELDS
SIM = orderfields(SIM);
end

%% PRE-SIMULATION TOOLS ///////////////////////////////////////////////////
% PREPARE THE PATCH FILES FOR VISUALS
function [SIM,STLimports]  = ConfigureVisuals(SIM)
% This function prepares object STL files for presentation. 
% ASSUMPTIONS:
% - The object stl file has the same name as the object being simulated.
% - The STL is correctly rotated to match the axes of XYZ of the agents
%   local frame.

% ABSOLUTE PATH TO THE OBJECT STL LOCATION
absPath = strcat(pwd,'\objects\');

% GENERATE FOR ALL OBJECTS
STLimports = 0;
for objectNumber = 1:SIM.totalObjects
    % GET THE AGENT SUPERCLASSES
    aliases = superclasses(SIM.OBJECTS(objectNumber).class);
    aliases = vertcat(SIM.OBJECTS(objectNumber).class,aliases);
    % MOVE THROUGH THE CLASS HEIRARCHY
    for aliasNo = 1:numel(aliases)
        % BUILD THE STL ASSOCIATION
        filename = strcat(absPath,char(aliases(aliasNo)),'.stl'); % Build the associated STL path
        % GET THE OBJECT PATCH FROM STL FILE (in the object directory) 
        [objectPatch,getFlag] = OMAS_visualTools.getSTLPatchObject(filename);
        if getFlag
           STLimports = STLimports + 1;
           break 
        end
    end
    
    % IF PATCH RETURNED, PROCESS FOR SIMULATION
    if isstruct(objectPatch)
        % GET MAXIMAL DIMENSION
        dimMaximal = max(max(abs(objectPatch.vertices),[],2),[],1);
        % CREATE DIMENSIONS SPECIFIED BY VIRTUAL.radius
        objectPatch.vertices = (objectPatch.vertices/dimMaximal)*SIM.OBJECTS(objectNumber).radius;
    end
    % RETURN PATCH-POPULATED META OBJECT
    SIM.OBJECTS(objectNumber).geometry = objectPatch;
end
end
% CONFIGURE THE PARALLEL/MULTITHREADED WORKERS
function [poolObject] = ConfigureThreadpool()
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
function [outputPath,fileString] = ConfigureOutputDirectory(MCenableFlag,absolutePath)
% INPUTS:
% workdir      - The current working directory string
% relOutputDir - The relative output directory 
% OUTPUTS:
% outpath - The current sessions output path
% filestr - The name of the session output file

% CONFIRM PATH CONVENTION
if ~strcmp(absolutePath(end),'\')
    absolutePath = strcat(absolutePath,'\');
end
absolutePath = OMAS_system.GetOSPathString(absolutePath);                  % Align path convention with OS

% Situations can occur where two simulation run at the same time, therefore
% we must handle the failed generation of the output directory.

% BUILD DESIRED OUTPUT PATH
fileString = [datestr(datetime('now'),'yyyy-mm-dd@HH-MM-SS'),'-session-data']; % Record current time                     % Build filestring

% DETERMINE IF MONTECARLO BEHAVIOUR IS REQUIRED
if MCenableFlag                             
    flag = 0;
    while flag ~= 1
        % RANDOMISE THE FILE NAME
        fileString = sprintf('[%s] cycle_data',sprintf('%d',randi([0,9],1,10)));
        flag = mkdir(absolutePath,fileString);                             % Create output directory
    end
    % Assign output path
    outputPath = OMAS_system.GetOSPathString([absolutePath,fileString,'\']);
    return
end
% Check the path is aligned with the host OS
outputPath = OMAS_system.GetOSPathString([absolutePath,fileString,'\']);   % Define the output directory
% Attempt to create new directory if required
if 7 ~= exist(outputPath,'file')
    success = mkdir(absolutePath,fileString);                              % Attempt to create the directory
    if ~success
        error('Output directory creation failed: %s',[absolutePath,fileString]);
    end
end
end