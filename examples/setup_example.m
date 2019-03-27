%% EXAMPLE SETUP (setup_example.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is designed as an example of how to simulate a multiple agent
% collision scenario using the Open Multi-Agent Simulation (OpenMAS) suite.

% For descriptions of how to create objects see:
% /objects/readme-objects.txt

% For a general overview see:
% /readme.txt

% Author: James A. Douthwaite 23/11/2018

clear all; close all;

%% ADD ADDITIONAL THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('scenarios');   

%% DEFINE SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%');              % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\OpenMAS_data');   % Set it as output path
sim_vebosity   = 1;
sim_warningDistance = 2;
sim_maxDuration = 10; 
sim_timeStep   =  0.1;           
% sim_figureSet = {'collisions'};         % SEE 'OMAS_figureGenerator.m' for options
sim_figureSet = {'events','plan','gif'};

%% SCENARIO PARAMETERS 
sim_agentNumber = 2;                   % TOTAL NUMBER OF AGENTS
sim_agentOrbit  = 5;                   
sim_waypointOrbit = 10;
sim_offsetAngle = 0;
sim_agentVelocity = 0;
sim_obstacleNumber = 4;
sim_noiseSigma = 0.2;
sim_plotScenario = logical(true);

%% GENERATE SCENARIO/OBSTACLE SET
fprintf('|| Assigning agent definitions:\n');
for index = 1:sim_agentNumber
%     agentIndex{index} = agent();
%     agentIndex{index} = agent_2D();

%     agentIndex{index} = agent_test();
    agentIndex{index} = agent_example();

% 3D COLLISION AVOIDANCE
%     agentIndex{index} = agent_VO();
%     agentIndex{index} = agent_RVO();
%     agentIndex{index} = agent_HRVO();
%     agentIndex{index} = agent_vectorSharing();

% 2D COLLISION AVOIDANCE
%     agentIndex{index} = agent_2D_VO();
%     agentIndex{index} = agent_2D_RVO();
%     agentIndex{index} = agent_2D_HRVO();
%     agentIndex{index} = agent_2D_vectorSharing();

% DYNAMICAL EXAMPLES
%     agentIndex{index} = ARdrone_LQR();
%     agentIndex{index} = quadcopter();
end

%% PLACE AGENT OBJECTS IN PRE-DEFINED SCENARIO
% OBSTACLE TESTS
% [ objectIndex ] = GetScenario_fourCuboidObstacles('agents',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_obstacleTrack('agents',agentIndex,'plot',sim_plotScenario);

% % AGENT TESTS
% [ objectIndex ] = GetScenario_twoLines('agents',agentIndex,'agentVelocity',sim_agentVelocity,'padding',2,'plot',sim_plotScenario);
[ objectIndex ] = GetScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',sim_offsetAngle,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = GetScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = GetScenario_concentricSphere('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_concentricAngle('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'angle',pi/4,'plot',sim_plotScenario);

% RANDOM TESTS
% [ objectIndex ] = GetScenario_random('objects',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_randomNormal('objects',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_randomUniform('objects',agentIndex,'plot',sim_plotScenario);

% WAYPOINT TESTS                            
% [ objectIndex ] = GetScenario_waypointCurve('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);
% [ objectIndex ] = GetScenario_waypoint_90Degrees('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);

% RL TESTS
% [ objectIndex ] = GetScenario_earthOrbit('plot',sim_plotScenario);

%% %%%%%% INITIALISE THE SIMULATION WITH THE OBJECT INDEX %%%%%%%%%%%%%%%%%
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,... 
                                   'dt',sim_timeStep,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
                           'outputPath',sim_outputPath);
clearvars -except DATA META
% /////////////////////////////////////////////////////////////////////////
% Upon completion of the simulation all the trajectory, figures and
% statistics will be pushed to the output directory, or held in 'DATA'. 

% IMPORT DATA (from 'outputPath')
load(strcat(META.outputPath,'META.mat'));
load(strcat(META.outputPath,'EVENTS.mat'));
fprintf('Complete...\n');