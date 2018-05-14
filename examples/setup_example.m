%% EXAMPLE SETUP (setup_example.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is designed as an example of how to simulate a multiple agent
% collision scenario using the Multi-Agent System Simulation (MASS) tool set.

% For descriptions of how to create objects see:
% /objects/readme-objects.txt

% For a general overview see:
% /readme.txt

% Author: James A. Douthwaite 16/05/2016

clear all; close all;

%% ADD THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('scenarios');   

%% 1. GENERATE SCENARIO/OBSTACLE SET
fprintf('|| Assigning agent definitions:\n');
agentNumber = 5;                                                           % The number of instances of the given test objects
for index = 1:agentNumber
%     agentIndex{index} = agent_test();
%     agentIndex{index} = agent_example();
    agentIndex{index} = agent_VO();
%     agentIndex{index} = agent_RVO();
%     agentIndex{index} = agent_HRVO();                                    % Try out some of the working instances
%     agentIndex{index} = agent_vectorSharing();                           % Define a cell array of test objects (full list available in /objects)
end

%% 2. GET A PRE-DEFINED CONCENTRIC SCENARIO
objectIndex = getScenario_concentric('agents',agentIndex,...               % The list of agent objects (agent_example,agent_example)
                                     'radius',50,...                       % The radius of the agent distribution
                                     'velocities',5,...                    % Define the initial forward velocity of the agents
                                     'plot',0);                            % Request the scenario be plotted

% BUILD GLOBAL OBJECT SET
clearvars -except objectIndex                                              % Clean up

%% 3. DEFINE THE SIMULATION PARAMETERS
% A full list of the available figures can be found in:
% /environment/simulation_figureIndex.m
% figureSet = {'all'};
figureSet = {'eventoverview','fig'}; 

%%\\\\\\ INITIALISE THE SIMULATION WITH THE OBJECT INDEX \\\\\\\\\\\\\\\\\\
% A full list of the available simulation options can be found in:
% \environment\simulation_initialise.m within 'validateSimulationInputs()'

[DATA,META] = simulation_initialise('objects',objectIndex,...              % Pass the list of initialised object/agents
                                    'simTime',30,...                       % Pass the max simulation time
                                         'dt',0.05,...                     % Pass the simulation step size
                                    'figures',figureSet);                  % Pass the requested figure list

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% clearvars -except DATA META
load(strcat(META.outputPath,'META.mat'));                                  % Load the META file that was saved
load(strcat(META.outputPath,'EVENTS.mat'));                                % Load the EVENTS file that was saved
fprintf('..Complete.\n');