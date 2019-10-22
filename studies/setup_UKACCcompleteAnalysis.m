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
monteCarlo_algs = {'agent_2D_VO',...
                   'agent_2D_RVO',...
                   'agent_2D_HRVO',...
                   'agent_2D_RVO2'};
% monteCarlo_algs = {'agent_2D_RVO2'};
monteCarlo_agentN = [2,5,10,20];
monteCarlo_noiseSigma = 0.2;
monteCarlo_outputDir = strcat(pwd,'\data\UKACC_realisticSensing');
monteCarlo_cycles = 1000;   % [FIXED]
monteCarlo_parallel = 1;
monteCarlo_offOnComplete = 0;

% FOR EACH NUMBER OF AGENTS IN THE GIVEN SCENARIO
for n = 1:numel(monteCarlo_agentN)                                         % Differing agent numbers 
    % FOR EACH ALGORITHM TO BE EVALUATED
    for index = 1:numel(monteCarlo_algs)                                   % Differing agent algorithms
        % CREATE THE AGENT SET
        agentVector = cell(1,monteCarlo_agentN(n));                        % Agent container
        for i = 1:monteCarlo_agentN(n)      
            agentVector{1,i} = eval(monteCarlo_algs{index});               % Assign agents
        end
        
        % ///////////////////// SCENARIO CONFIGURATION ////////////////////
        initialisedObjects = getScenario_UKACC(...                         % Initialise the objects in the UKACC scenario
            'agents',agentVector,...
            'agentOrbit',sim_agentOrbit,...
            'agentVelocity',sim_agentVelocity,...
            'waypointOrbit',sim_waypointOrbit,...
            'offsetAngle',0);
        clear agentVector
        
        % UNIQUE OUTPUT FILE
        studyOutputDir = strcat(monteCarlo_outputDir,...
                                '\',monteCarlo_algs{index});
        
        % ///////////////// DEFINE THE MONTE-CARLO INSTANCE ///////////////
        fprintf('[SETUP]\t Defining Monte-Carlo Simulation Series...\n');
        [monteCarlo] = OMAS_monteCarlo(...
                        'objects',initialisedObjects,...
                        'duration',sim_maxDuration,...
                        'dt',sim_timeStep,...
                        'warningDistance',sim_warningDistance,...
                        'verbosity',sim_verbosity,...
                        'outputPath',studyOutputDir,...
                        'positionSigma',monteCarlo_noiseSigma,...
                        'threadPool',monteCarlo_parallel,...
                        'shutDownOnComplete',monteCarlo_offOnComplete,...
                        'cycles',monteCarlo_cycles);

        % ///////////////// EXECUTE CYCLE EVALUATION //////////////////////
        [summaryData,monteCarlo] = monteCarlo.evaluateCycles();            % Process the defined cycles
        clear summaryData monteCarlo
    end    
end               