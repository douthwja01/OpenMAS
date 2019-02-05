%% SOLAR SYSTEM EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This setup function configures an objectIndex that resembles our solar
% system in their orbital configuration. 

% ADD THE PROGRAM PATHS
clear all; close all; 
addpath('environment');
addpath('objects');  
addpath('scenarios'); 

fprintf('[SETUP]\tAssembling the solar system example.\n');

%% SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%'); % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\OpenMAS_data');
sim_vebosity   = 1;
sim_warningDistance = 1;
sim_threadPool = 0;
sim_maxDuration = 92*60; 
sim_timeStep    =  60;      
sim_figureSet = {'fig','gif'}; 

%% SCENARIO PARAMETERS 
sim_plotScenario = 1;
sim_agentOrbit  = 3; 
sim_waypointOrbit = 10;
sim_offsetAngle = 0;
sim_agentVelocity = 0;
sim_noiseSigma = 0.0;

%% INITIALISE THE OBJECT SCENARIO
fprintf('[SETUP]\tAssigning object definitions:\n');

objectIndex = getScenario_earthOrbit('plot',sim_plotScenario,'scale',5E3);                 % Generate the earth


%% %%%%%% INITIALISE THE SIMULATION WITH THE OBJECT INDEX %%%%%%%%%%%%%%%%%
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,... 
                                   'dt',sim_timeStep,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
                           'threadPool',sim_threadPool,...
                           'outputPath',sim_outputPath);
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except DATA META
load(strcat(META.outputPath,'META.mat'));
load(strcat(META.outputPath,'EVENTS.mat'));