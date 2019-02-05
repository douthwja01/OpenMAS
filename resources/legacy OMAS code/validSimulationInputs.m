% ERROR CHECKING //////////////////////////////////////////////
% ADD SOME DEFAULT PARAMETERS TO THE 'defaultConfig'
if ~isfield(defaultConfig,'file') || isempty(defaultConfig.file)
    defaultConfig.file = 'scenario.mat';
end
% CHECK AGENTS HAVE BEEN PROVIDED
if ~isfield(defaultConfig,'agents') || isempty(defaultConfig.agents)
    error("[SCENARIO]\t You must specify and agent cell set using the 'agents' input parameter");
end


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
                 'verboseMode',2,...                                       % Default to displaying debug information
               'visualisation',boolean(0));                                % Output visualisation
             
% DEFAULT TIMING PARAMETERS
SIM.TIME = struct('simTime',10,'dt',0.1,'startTime',0);                    % Assign default timing parameters        
SIM.TIME.maxRunTime = SIM.TIME.simTime;                                    % Determine maximum run time

%% BEGIN PARSING INPUT VARIABLES  /////////////////////////////////////////         
fprintf('[SETUP]\tConfirming input variables.\n');

%% TIMING PARAMETERS //////////////////////////////////////////////////////
tmp = strncmpi(inputSet,'SIMTIME',7)|strncmpi(inputSet,'RUNTIME',7); 
if any(tmp)
    SIM.TIME.simTime = inputSet{find(tmp) + 1};
    assert(isnumeric(SIM.TIME.simTime) == 1,'The simulation total time "simTime" must be numeric.');
end
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
    value = inputSet{find(tmp) + 1};
    assert(isnumeric(value) == 1,'Verbosity must be specified as a number between 0 and 2.');
    SIM.verboseMode = value; % Assume a value higher than one is a request
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