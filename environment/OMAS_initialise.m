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

% ///////////////// ASSIGN DEFAULT SIMULATION PARAMETERS //////////////////
% The user must be able to run a simulation for an intended amount of time,
% or run with a maximal time for all goals to be complete. 

% /////////////////////// DEFAULT TIMING STRUCTURE ////////////////////////
TIME = struct('duration',1E3,...                                           % Assign default simulation duration (for output matrix scaling)
                    'dt',0.5,...
             'startTime',0);
         
% //////////////// DEFAULT GENERAL CONFIGURATION STRUCTURE ////////////////
defaultConfig = struct('figures',{{'events','fig'}},...
               'warningDistance',10,...
            'conditionTolerance',1E-3,...
            'visabilityModifier',1,...
                    'outputPath',strcat(pwd,'\data'),... 
                    'systemFile','temp.mat',...
                         'phase','SETUP',...                               % Initialise in setup mode
                    'threadPool',logical(false),...
                'monteCarloMode',logical(false),...
                     'verbosity',2,...
                           'gui',logical(true),...
                       'objects',[],...
                          'TIME',TIME);
                 
% ///////////////////// PARSE USER INPUT PARAMETERS ///////////////////////
fprintf('[SETUP]\tConfirming input variables.\n');
[SIM] = parseConfigurationVector(defaultConfig,varargin);
clear TIME defaultConfig

% /////////////// SOME PARAMETER DEDUCTION/INPUT HANDLING /////////////////
assert(~isempty(SIM.objects),'Please provide a agent/object object vector.');
% TRIM ANY EMPTY OBJECT CONTAINERS
for ind = 1:numel(SIM.objects)
   if ~isobject(SIM.objects{ind})
      SIM.objects(ind) = {''};                                             % Remove invalid classes from the container
   end
end
SIM.objects = SIM.objects(~cellfun('isempty', SIM.objects)); 

% DISABLE AGENT-WORKER ALLOCATION IN MONTE-CARLO MODE
if SIM.monteCarloMode 
    SIM.threadPool = logical(false);
end

% SIMULATION PHASE INDICATOR
fprintf('[%s] --> Input parameters defined.\n',SIM.phase);

%% SIMULATION CONFIGURATION PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[%s]\tINITIALISING SIMULATION...\n\n',SIM.phase);
close all; 

% GET THE SIMULATIONS SUBDIRECTORIES (if not already added to path)
configureDependencies(pwd,{'environment','toolboxes','events','objects','scenarios'});

% MAKE ANY ADDITIONAL PREPERATIONS TO THE META DATA ///////////////////////
fprintf('[%s]\tOMAS CONFIGURATION:\n',SIM.phase);
clear META
[META,objectIndex] = defineSimulationParameters(SIM);
clearvars SIM;
fprintf('[%s]\tObjects parameterized.\n',META.phase);

% INITIALISE THE GUI (IF EXISTS AND IS REQUESTED) /////////////////////////
fprintf('[%s]\tGRAPHICS CONFIGURATION:\n',META.phase);
fprintf('[%s]\tImporting object STLs.\n',META.phase);
if META.gui && ~META.monteCarloMode
   [META,STLimports] = configureVisuals(META);
   fprintf('[%s]\t...%s files imported.\n',META.phase,num2str(STLimports));
   clear STLimports;
else
   fprintf('[%s]\t... STL importation supressed.\n',META.phase); 
end
% CONFIGURE THE OUTPUT SETTINGS ///////////////////////////////////////////
fprintf('[%s]\tOUTPUT CONFIGURATION:\n',META.phase);
% DEFINE THE RELATIVE OUTPUT PATH
[META.outputPath,META.outputFile] = configureOutputDirectory(META.monteCarloMode,META.outputPath);  % Get current working working directory variables
fprintf('[%s]\tSession output directory created: \n[%s]\t"%s" \n',META.phase,META.phase,META.outputFile);
fprintf('[%s]\tOutput path: \n[%s]\t"%s"\n',META.phase,META.phase,META.outputPath); 

% SAVE THE META DATA OBJECT
save([META.outputPath,'META.mat'],'META'); % Record the simulation parameters

% PUSH FIGURE REQUEST TO FILE
figureList = META.figures;
save([META.outputPath,META.systemFile],'figureList');
META = rmfield(META,'figures');            % Remove the figure list from the META structure
clear figureList

% DISPLAY SIMULATION CONFIGURATION 
fprintf('[%s]\tCONFIGURATION SUMMARY:\n',META.phase);
fprintf('[%s]\tObjects: %s\tAgents: %s\tObstacles: %s\tWaypoints: %s\n',...
        META.phase,num2str(META.totalObjects),num2str(META.totalAgents),...
        num2str(META.totalObstacles),num2str(META.totalWaypoints));                      % Object summary
fprintf('[%s]\tDuration: %ss\tSampling: %ss\tSteps: %s\n',META.phase,...
        num2str(META.TIME.duration),num2str(META.TIME.dt),num2str(META.TIME.numSteps)); % Timing summary
   
%% //////// BEGIN RUNNING THROUGH THE SIMULATION TIME CYCLE ///////////////
fprintf('[%s]\tMoving to simulation...\n',META.phase);
[DATA,META,EVENTS,objectIndex] = OMAS_mainCycle(META,objectIndex);

%% ////////////////////////////////////////////////////////////////////////

% BEGIN RUNNING EXTERNAL DATA ANALYSIS
% try
    fprintf('[%s]\tMoving to post-simulation analysis...\n[%s]\n',META.phase,META.phase);
    [DATA] = OMAS_analysis(META,objectIndex,EVENTS,DATA);                  % Jump to external analysis program
% catch analysisError
%     warning('[ERROR] A problem occurred interpreting the output data');
%     warning(analysisError.message);
% end
fprintf('\n[%s]\tExiting...\n\n',META.phase);
end

%% SIMULATION SETUP FUNCTIONS /////////////////////////////////////////////
% PREPARE THE META STRUCTURE AND INITIALISE OBJECT INDEX
function [SIM,objectIndex] = defineSimulationParameters(SIM)
% INPUTS:
% SIM         - The simulation inputs
% objectIndex - A vector of object classes
% OUTPUTS:
% SIM         - The prepared simulations parameters

% GENERAL PARAMETER INITIALISATION
% DETERMINE IF A PARALLEL POOL IS NEEDED
if SIM.threadPool 
    [SIM.threadPool] = configureThreadpool();
end

% ////////////// SIMULATION PARAMETER INITIALISATION //////////////////////
% DEFINE TIME PARAMETERS
% SIM.TIME.maxTime    = SIM.TIME.startTime + SIM.TIME.simTime;               % Define the end time
SIM.TIME.endTime    = SIM.TIME.duration + SIM.TIME.startTime;              
SIM.TIME.timeVector = SIM.TIME.startTime:SIM.TIME.dt:(SIM.TIME.endTime);   % The progressive time vector
SIM.TIME.numSteps   = size(SIM.TIME.timeVector,2);                         % Number of resulting steps
SIM.TIME.frequency  = 1/SIM.TIME.dt;                                       % Define the simulation frequency 

% INSTANTEOUS PARAMETERS
SIM.TIME.currentTime = 0;
SIM.TIME.currentStep = 1;

% DEFINE THE AVAILABLE STATUS TYPES FROM THE 'eventType' ENUMERATIONS
% [~,SIM.eventNames] = enumeration('eventType');
% SIM.eventNames = SIM.eventNames(2:end);
% SIM.eventTypes = length(SIM.eventNames);                                   % Size the number of available sim events (neglect base definition)
[eventEnums,~] = enumeration('eventType');

% SEPERATE THE OBJECT VECTOR
objectIndex = SIM.objects; 
SIM = rmfield(SIM,'objects');                                              % Remove cell array to reduce data overhead

% DEFINE THE NUMBER OF AGENTS
SIM.totalObjects = numel(objectIndex);                                     % Define the object number

%PREPARE THE META.OBJECT GLOBAL STATUS VECTOR
for i = 1:SIM.totalObjects
    % THE OBJECTS ARE ASSIGNED THEIR GLOBAL PARAMETERS IN THE GLOBAL XYZ 
    % (MATLAB) AXIS SYSTEM;
    % .position     - Global position in XYZ coordinates 
    % .velocity     - Global velocity in XYZ coordinates
    % .quaternion   - Quaternion describing the XYZ rotation
    
    % THE SIMULATION OPERATES IN THE ENU FROM IDENTICAL TO THE MATLAB FRAME
    
    % GLOBAL PARAMETERS ARE DECLARED IN THE MATLAB XYZ FRAME
    globalENUPosition = objectIndex{i}.VIRTUAL.globalPosition;              % Get ENU global position
    globalENUVelocity = objectIndex{i}.VIRTUAL.globalVelocity;              % Get ENU global velocity
    staticBodyToRotatedGlobalQuaternion = objectIndex{i}.VIRTUAL.quaternion;% REDEFINE INITIAL QUATERNION AS STATIC BODY TO ROTATED EARTH QUATERNION
    
    % CONFIRM XYZ INPUTS
    assert(numel(globalENUPosition) == 3,'Global position must be given as a vector [3x1].');
    assert(numel(globalENUVelocity) == 3,'Global velocity must be given as a vector [3x1].');
    assert(numel(staticBodyToRotatedGlobalQuaternion) == 4,'Global attitude must be give as a quaternion vector [4x1].');      
    % GET THE INITIAL ROTATIONS (XYZ -> XYZ)
    [R_BG,R_GB] = OMAS_axisTools.quaternionToRotationMatrix(staticBodyToRotatedGlobalQuaternion);
    % CALCULATE THE EULER HEADING
    ZYXrotations = rotm2eul(R_GB,'ZYX');
    ENUrotations = [0,0,1;0,1,0;1,0,0]*ZYXrotations';                      % Get the equivalent rotations around the XYZ axes
        
    % /////// UNTIL THIS POINT, GLOBAL PARAMETERS ARE GENERALISED /////////
    % The problem we find is that we don't know how the objects states 
    % relate to the allocated global position, velocity and quaternion.   
    % Instead we call the objects method for initialising its own state.
    objectIndex{i} = objectIndex{i}.initialise_localState(R_GB*globalENUVelocity,ENUrotations);

    % ///////////////// INDIVIDUAL OBJECT ASSESSEMENT /////////////////////
    % ASSIGN THE INITIAL ROTATION MATRIX
    if ~isfield(objectIndex{i}.VIRTUAL,'R')
        objectIndex{i}.VIRTUAL.R = R_BG;
    end
    % IF THERE IS NO DETECTION RANGE GIVEN (object/blind agent)
    if ~isfield(objectIndex{i}.VIRTUAL,'detectionRange')
        objectIndex{i}.VIRTUAL.detectionRange = 0;
    end
    
    % ///// CHECK FOR FORBIDDEN CHARACTERS IN OBJECT NAME ASSIGNMENTS /////
    % This prevents issues from occurring with filenames and figures.
    assertionString = sprintf('Forbidden character (:*?"<>|) in object identifier: "%s"',objectIndex{i}.name);
    assert(~any(contains(objectIndex{i}.name,{'\','/',':','*','?','"','<','>','|'},'IgnoreCase',true)),assertionString);
    
    % OBJECTS is a vector of structures containing the object META data 
    SIM.OBJECTS(i) = struct('objectID',objectIndex{i}.objectID,...                  % ID reference number
                                'name',objectIndex{i}.name,...                      % Name reference
                               'class',class(objectIndex{i}),...                    % Class reference
                                'type',objectIndex{i}.VIRTUAL.type,...              % Object type
                              'colour',objectIndex{i}.VIRTUAL.colour,...            % Object virtual colour
                              'symbol',objectIndex{i}.VIRTUAL.symbol,...            % Representative symbol (from type)
                               'patch',[],...                                       % The objects patch verticies and faces
                              'radius',objectIndex{i}.VIRTUAL.radius,...            % The objects critical radius
                      'detectionRange',objectIndex{i}.VIRTUAL.detectionRange,...    % The simulations detection horizon for the aircraft
                          'idleStatus',objectIndex{i}.VIRTUAL.idleStatus,...        % Agent completed Task Flag
                         'globalState',[globalENUPosition;globalENUVelocity;staticBodyToRotatedGlobalQuaternion],...
                                'R_GB',R_GB,...                                     % Current Global-body rotation matrix
                                'R_BG',R_BG,...                                     % Current Body-global rotation matrix
                   'relativePositions',zeros(SIM.totalObjects,3),...              % Relative distances between objects [dx;dy;dz]*objectCount
                 'euclideanSeparation',zeros(SIM.totalObjects,1),...                % Euclidean scalar distances between objects (excluding itself)
                        'objectStatus',zeros(SIM.totalObjects,(numel(eventEnums)-1)/2));
    % DEFINE THE OBJECT-STATUS MATRIX
%     SIM.OBJECTS(i).objectStatus = zeros((SIM.totalObjects),(numel(eventEnums/2)));  % A matrix recording all interactions against the object
                    
    % ASSIGN NaN TO LOCATIONS OF SELF-REFERENCE            
    SIM.OBJECTS(i).relativePositions((objectIndex{i}.objectID),:) = NaN;
    SIM.OBJECTS(i).euclideanSeparation((objectIndex{i}.objectID),1) = NaN;
    SIM.OBJECTS(i).objectStatus((objectIndex{i}.objectID),:) = NaN;
  
    % CHECK FOR SIMULATION/SAMPLING FREQUENCY ERROR
    % The user specifies an agent control frequency, which must be lower
    % than the minimum simulation sample frequency (1/dt)
    if isprop(objectIndex{i},'SENSORS') && isfield(objectIndex{i}.SENSORS,'sampleFrequency')
        if isinf(objectIndex{i}.SENSORS.sampleFrequency)
        	objectIndex{i}.SENSORS.sampleFrequency = (1/SIM.TIME.dt);
        elseif objectIndex{i}.SENSORS.sampleFrequency  > (1/SIM.TIME.dt)
            error('Agent %s sample frequency (%sHz) cannot be greater than the simulation frequency (%sHz).',...
                    objectIndex{i}.name,num2str(objectIndex{i}.sampleFrequency),num2str(1/SIM.TIME.dt));
        end
    end
end

% GET THE TOTALS FOR EACH OBJECT TYPE FOR REFERENCE
[objEnums,objNames] = enumeration('OMAS_objectType');                      % Define the number of other object types
typeIndex = [SIM.OBJECTS.type];
for i = 1:numel(objEnums)
    typeLabel = strcat('total',char(objNames(i)),'s');
    typeLabel(6) = upper(typeLabel(6));
    SIM.(typeLabel) = sum(typeIndex == objEnums(i));
end
% GET THE OBJECT ID REFERENCE VECTOR
SIM.globalIDvector = [SIM.OBJECTS.objectID];
% ORGANISE THE FIELDS
SIM = orderfields(SIM);
end

% PREPARE THE PATCH FILES FOR VISUALS
function [SIM,STLimports] = configureVisuals(SIM)
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
    SIM.OBJECTS(objectNumber).patch = objectPatch;
end
end

%% PRE-SIMULATION TOOLS ///////////////////////////////////////////////////
% PARSE A PARAMETER VECTOR AGAINST A CONFIGURATION STRUCTURE
function [config] = parseConfigurationVector(defaultConfig,parameterVector)
% This function is designed to parse a generic set of user
% inputs and allow them to be compared to a default input
% structure. The function searches for associated properties and
% re-allocates the named variable with the associated value.

% MOVE THROUGH THE PARAMETER PAIRS ('agents',agentIndex)
for parameterIndex = 1:numel(parameterVector)
    % FOR EACH USER-DEFINED PARAMETER
    givenParameter = parameterVector{parameterIndex};
    if ~ischar(givenParameter)
        continue
    end
    assumedValue = parameterVector{parameterIndex + 1}; % The associated value
    
    % IF THE OBJECT HAS A PROPERTY BY THAT NAME
    if isfield(defaultConfig,givenParameter)
        defaultConfig.(givenParameter) = assumedValue;   % Make a substitution
    end
    % IF THE CONFIG HAS A SUBSTRUCTURE WITH THAT PROPERTY
    if isfield(defaultConfig.TIME,givenParameter)
        defaultConfig.TIME.(givenParameter) = assumedValue;
    end
end
% REASSIGN FOR CLARITY
config = defaultConfig;
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
% GET THE SIMULATION SUBDIRECTORY SET
function configureDependencies(wrkDir,directorySet)
    % Check the directory set is valid
    if ~iscell(directorySet)
        warning('Simulation subpath vector is invalid.');
        return
    end   
    
    % Move through the sub-directory list
%     for i = 1:size(directorySet,2)
%         pathName = strcat(wrkDir,'\',directorySet{i});
%         try
%             addpath(pathName);
%         catch simDirectoryError
%             warning('Unable to find critical sub-directories, have you changed directories?');
%             error(simDirectoryError.message);
%         end
%     end

pathCell = regexp(path, pathsep, 'split');

if ispc  % Windows is not case-sensitive
    for i = 1:size(directorySet,2)
        pathName = strcat(wrkDir,'\',directorySet{i});
        if ~any(strcmpi(pathName, pathCell))
            addpath(pathName);
        end
    end
else
    for i = 1:size(directorySet,2)
        pathName = strcat(wrkDir,'\',directorySet{i});
        if ~any(strcmp(pathName, pathCell))
            addpath(pathName);
        end
    end
end
end
% CONFIGURE THE OUTPUT DIRECTORY
function [outputPath,fileString] = configureOutputDirectory(MCenableFlag,absolutePath)
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

% Situations can occur where two simulation run at the same time, therefore
% we must handle the failed generation of the output directory.

% BUILD DESIRED OUTPUT PATH
fileString = strcat('sessiondata',datestr(datetime('now'),' yyyy-mm-dd @ HH-MM-SS')); % Record current time                     % Build filestring
outputPath = strcat(absolutePath,fileString,'\');                          % Define the output directory

% DETERMINE IF MONTECARLO BEHAVIOUR IS REQUIRED
if MCenableFlag                             % [TO BE REVISED]
    flag = 0;
    while flag ~= 1
        % RANDOMISE THE FILE NAME
        fileString = strcat('sessiondata [',sprintf('%d',randi([0,9],1,10)),']');
        outputPath = strcat(absolutePath,fileString,'\');                  % Define the output directory
        % CREATE OUTPUT DIRECTORY
        flag = mkdir(absolutePath,fileString);
    end
    return
end

if 7 ~= exist(outputPath,'file')
    assert(mkdir(absolutePath,fileString) == 1,'Output directory creation failed.');
    return
else
    % Assume we intend to output the data to the exising directory
end
end