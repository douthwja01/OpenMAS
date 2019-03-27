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

fprintf('[SETUP]\tInitialising Monte-Carlo example.\n');
% /////////////////// AGENT & SCENARIO PARAMETERS /////////////////////////
sim_agentOrbit  = 30;       
sim_agentVelocity = 0;      
sim_waypointOrbit = 35;     

% //////////////////////// OMAS CONFIGURATION /////////////////////////////
sim_warningDistance = 10;   
sim_maxDuration = 120;      
sim_timeStep = 0.25;     	
sim_verbosity = 0;          

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
monteCarlo_algorithms = {'agent_2D_VO',...
                         'agent_2D_RVO',...
                         'agent_2D_HRVO'};
                         %'agent_2D_RVO2'};                                % The types of agents simulated 
monteCarlo_agentPopulations = [2,4]; 

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
[~, userdir] = system('echo %USERPROFILE%');    % Get desktop path
monteCarlo_outputPath = strcat(userdir,'\desktop\US18_data');    
monteCarlo_positionSigma = 0.2;
monteCarlo_velocitySigma = 0;
monteCarlo_cycles = 5;                        % [FIXED]
monteCarlo_parallel = 0;        
monteCarlo_offOnComplete = 0;

% ///////////////////// DEFINE THE MC OBJECT MATRIX ///////////////////////
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
        % Concentric scenario
%         initialisedObjects = GetScenario_concentric('agents',agentVector,...
%                                                     'agentOrbit',sim_agentOrbit,...
%                                                     'agentVelocity',sim_agentVelocity,...
%                                                     'waypointOrbit',sim_waypointOrbit);
        % Opposing lines scenarios                                             
        initialisedObjects = GetScenario_twoLines('agents',agentVector,...
                                                  'agentVelocity',sim_agentVelocity);                                             
        % /////////////////////////////////////////////////////////////////
        sessionSample = studyVector;
        sessionSample(1:numel(initialisedObjects)) = initialisedObjects;
        studyObjects = vertcat(studyObjects,sessionSample); 
    end
end

%% ///////////// BEGIN SIMULATING THE MONTE-CARLO INSTANCES ///////////////
fprintf('[SETUP]\t Defining Monte-Carlo Simulation Series...\n');
[monteCarlo] = OMAS_monteCarlo('objects',studyObjects,...
    'duration',sim_maxDuration,...
    'dt',sim_timeStep,...
    'warningDistance',sim_warningDistance,...
    'verbosity',sim_verbosity,...
    'positionSigma',monteCarlo_positionSigma,...
    'velocitySigma',monteCarlo_velocitySigma,...
    'threadPool',monteCarlo_parallel,...
    'shutDownOnComplete',monteCarlo_offOnComplete,...
    'cycles',monteCarlo_cycles,...
    'outputPath',monteCarlo_outputPath);

% EVALUATE ALL OF THE PROPOSED CONDITIONS
monteCarlo.evaluateAllCycles();  