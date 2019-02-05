%% EXAMPLE MONTE-CARLO SETUP CONFIGURATION FILE (monteCarlo_example.m) %%%%
% This function is intended to demonstrate how the monte-carlo aspect of
% OMAS may be used to run numerous cycles of a series of different
% algorithms, agents and scenarios.

% Author: James A. Douthwaite 30/04/18

% Notes*
% - The monte-carlo object accepts a cell array of objects [studies,objects]
%   (where 'studies' are seen as independant sets of cycles/analytics).
% - It is currently more efficient to hand smaller sets of object arrays 
%   per Monte-Carlo instance rather than run many studies using one.
% - The object array can be defined in same manner as OMAS, using any 
%   scenario or object set, the Monte-Carlo instance will then perturb the
%   global positions of all objects before the cycle begins.

% ADD THE PROGRAM PATHS
clear all; close all;
addpath('environment');
addpath('objects');
addpath('scenarios');                                                      % Required system/toolbox paths

% /////////////////// AGENT & SCENARIO PARAMETERS /////////////////////////
fprintf('[SETUP]\tInitialising Monte-Carlo example.\n');
scenario_agentOrbit  = 30;       
scenario_agentVelocity = 0;      
scenario_waypointOrbit = 35;                                               % Parameters used to define the scenario

% //////////////////////// OMAS CONFIGURATION /////////////////////////////
sim_warningDistance = 10; 
sim_maxDuration = 120;    
sim_timeStep = 0.25;      
sim_verbosity = 0;                                                         % Parameters defining the OMAS configuration within a give cycle. 

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
monteCarlo_algs = {'agent_2D_VO','agent_2D_RVO','agent_2D_HRVO','agent_2D_RVO2'}; % The selected agents within the '/objects' folder
monteCarlo_agentN = [2,4,6,8];                                             % Run cycles for 2, 4, 6, 8 agents in sequence
monteCarlo_noiseSigma = 0.2;                                               % The position purturbation of the objects
monteCarlo_outputDir = strcat(pwd,'\data');                                % The directory to output the MC and OMAS data
monteCarlo_cycles = 10;                                                    % The number of cycles per study
monteCarlo_parallel = 1;                                                   % Use the Parallel Toolbox to handle the cycles
monteCarlo_offOnComplete = 0;                                              % Turn off the PC on study completion.

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
            'agentOrbit',scenario_agentOrbit,...
            'agentVelocity',scenario_agentVelocity,...
            'waypointOrbit',scenario_waypointOrbit,...
            'offsetAngle',0);
        clear agentVector
        
        % UNIQUE OUTPUT FILE
        studyOutputDir = strcat(monteCarlo_outputDir,...                   % Output to path defined by the object name
                                '\',monteCarlo_algs{index});    
        
        % ///////////////// DEFINE THE MONTE-CARLO INSTANCE ///////////////
        fprintf('[SETUP]\t Defining Monte-Carlo Simulation Series...\n');   
        [monteCarlo] = OMAS_monteCarlo(...                                 % Create an instance of the Monte-Carlo object
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
        [summaryData,monteCarlo] = monteCarlo.evaluateCycles();            % Evaluate the scenarios for the given number of cycles
        clear summaryData monteCarlo                                       % Clear memory
    end    
end               