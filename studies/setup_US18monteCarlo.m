% Author: James A. Douthwaite 26/03/18

% ADD THE PROGRAM PATHS
clear all; close all;
addpath('environment');
addpath('objects');
addpath('scenarios');

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
   