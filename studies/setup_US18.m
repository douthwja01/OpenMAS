%% OpenMAS MONTE-CARLO SIMULATION SETUP (setup_US18.m) %%%%%%%%%%%%%%
% This function is designed to build a monte-carlo style analysis of a
% given scenario. 

% Author: James A. Douthwaite 26/03/18

% ADD THE PROGRAM PATHS
clear all; close all; 
addpath('environment');
addpath('objects');  
addpath('scenarios'); 
addpath('toolboxes/Intlab_V7.1');   

%If intlab needs to be reloaded
try 
    wrkDir = pwd;
    test = infsup(0,1);
    clearvars test 
    IntDir = strcat(pwd,'\Intlab_V7.1\startupJD.m');
catch 
    run('startup.m'); 
    cd(wrkDir)
end

%% //////////////////// AGENT & SCENARIO PARAMETERS ///////////////////////
fprintf('[SETUP]\tInitialising Monte-Carlo setup script.\n');

%% SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%'); % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\OpenMAS_data');
sim_warningDistance = 10;
sim_threadPool = 0;
sim_maxDuration = 30;
sim_timeStep = 0.25;    % RVO timestep (0.25s)
sim_vebosity = 1;
% sim_figureSet ={'all'};
sim_figureSet = {'plan','closest','inputs','gif'}; 
%{'fig','gif','inputs','separations'};

% SCENARIO PARAMETERS
sim_plotScenario = 1;
sim_agentNumber = 8;    % TOTAL NUMBER OF AGENTS
sim_agentOrbit  = 10;
sim_waypointOrbit = 10;
sim_agentVelocity = 0;
sim_noiseSigma = 0.2;

% ////////////////////////// INITIALISE AGENTS ////////////////////////////
fprintf('[SETUP]\tAssigning agent definitions:\n');
for index = 1:sim_agentNumber
    % BASIC CLASSES
%     agentIndex{index} = agent_VO();
%     agentIndex{index} = agent_RVO();
%     agentIndex{index} = agent_HRVO();
%     agentIndex{index} = agent_vectorSharing();
%     agentIndex{index} = agent_vectorSharing_interval();
  
%     agentIndex{index} = agent_2D_VO();
    agentIndex{index} = agent_2D_RVO();
%     agentIndex{index} = agent_2D_HRVO();
%     agentIndex{index} = agent_2D_RVO2(); 
%     agentIndex{index} = agent_2D_vectorSharing();
%     agentIndex{index} = agent_2D_vectorSharing_interval(); 
end

fprintf('[SETUP]\t Intialising the agents in the given scenario.\n');
% CALL THE FUNCTION THAT GENERATES THE UKACC SCENARIO
% [ objectIndex ] = getScenario_UKACC_concentric('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario);
[ objectIndex ] = getScenario_US_twoLines('agents',agentIndex,'agentVelocity',sim_agentVelocity,'noiseFactor',sim_noiseSigma,'plot',sim_plotScenario);

% ////////////// INITIALISE THE SIMULATION WITH THE OBJECT INDEX /////////////////
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,... 
                                   'dt',sim_timeStep,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
                           'threadPool',sim_threadPool,...
                           'outputPath',sim_outputPath);
                       
clearvars -except DATA META
load(strcat(META.outputPath,'META.mat'));
load(strcat(META.outputPath,'EVENTS.mat'));


% ===================== MONTE-CARLO SIMULATION ============================
run_monteCarlo = false;
if ~run_monteCarlo 
    return
end
    
fprintf('[SETUP]\tInitialising Monte-Carlo setup script.\n');
% /////////////////// AGENT & SCENARIO PARAMETERS /////////////////////////
sim_agentOrbit  = 30;       % [FIXED]
sim_agentVelocity = 0;      % [FIXED]
sim_waypointOrbit = 35;     % [FIXED]

% //////////////////////// OMAS CONFIGURATION /////////////////////////////
sim_warningDistance = 10;   % [FIXED]
sim_maxDuration = 120;      % [FIXED]
sim_timeStep = 0.25;     	% [FIXED] 
sim_verbosity = 0;          % [FIXED]

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
monteCarlo_algorithms = {'agent_2D_VO',...
                         'agent_2D_RVO',...
                         'agent_2D_HRVO'};
                         %'agent_2D_RVO2'};     % The types of agents simulated 
monteCarlo_agentPopulations = [2,5,10]; 

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
[~, userdir] = system('echo %USERPROFILE%');    % Get desktop path
monteCarlo_outputPath = strcat(userdir,'\desktop\US18_data');    
monteCarlo_noiseSigma = 0.2;
monteCarlo_cycles = 5;                        % [FIXED]
monteCarlo_parallel = 1;        
monteCarlo_offOnComplete = 0;

% /////////////////////////////////////////////////////////////////////////
studyVector = cell(1,2*max(monteCarlo_agentPopulations));
studyObjects = [];
for n = 1:numel(monteCarlo_agentPopulations)                               % Differing agent numbers 
    % FOR EACH ALGORITHM TO BE EVALUATED
    for index = 1:numel(monteCarlo_algorithms)                             % Differing agent algorithms
        % CREATE THE AGENT SET
        agentVector = cell(1,monteCarlo_agentPopulations(n));              % Agent container
        for i = 1:monteCarlo_agentPopulations(n)      
            agentVector{1,i} = eval(monteCarlo_algorithms{index});         % Assign agents
        end
        % ///////////////////// SCENARIO CONFIGURATION ////////////////////
        initialisedObjects = getScenario_UKACC('agents',agentVector,...
                                               'agentOrbit',sim_agentOrbit,...
                                               'agentVelocity',sim_agentVelocity,...
                                               'waypointOrbit',sim_waypointOrbit,...
                                               'offsetAngle',0);
        sessionSample = studyVector;
        sessionSample(1:numel(initialisedObjects)) = initialisedObjects;
        studyObjects = vertcat(studyObjects,sessionSample); 
    end
end

%% //////////////// CREATE THE MONTE-CARLO INSTANCE ///////////////////////
fprintf('[SETUP]\t Defining Monte-Carlo Simulation Series...\n');
[monteCarlo] = OMAS_monteCarlo('objects',studyObjects,...
    'duration',sim_maxDuration,...
    'dt',sim_timeStep,...
    'warningDistance',sim_warningDistance,...
    'verbosity',sim_verbosity,...
    'positionSigma',monteCarlo_noiseSigma,...
    'threadPool',monteCarlo_parallel,...
    'shutDownOnComplete',monteCarlo_offOnComplete,...
    'cycles',monteCarlo_cycles,...
    'outputPath',monteCarlo_outputPath);

% EVALUATE ALL OF THE PROPOSED CONDITIONS
monteCarlo.evaluateAllCycles();
