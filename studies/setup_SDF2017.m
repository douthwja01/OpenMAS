%% SENSOR DATA FUSION 2017 (setup_SDF2017.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the experimental data used in the SDF2017 paper.
% This configuration is not intended to be changed, only to run these
% specific tests to reproduce example figures

% Test data will be outputted to the three most recent entries in the
% 'data' directory.

clear all; close all;

% GET THE CURRENT DIRECTORY DETAILS
parentDir = 'studies';
cd ..
% ADD SIMULATION PATHS
addpath('environment');        
addpath('objects');  
addpath('scenarios'); 
addpath('Intlab_V7.1');   

% LOAD THE INTERVAL ELEMENTS
wrkDir = pwd;
try 
    test = infsup(0,1);
    clearvars test 
    IntDir = strcat(pwd,'\Intlab_V7.1\startupJD.m');
catch 
    run('startup.m'); 
    cd(wrkDir)
end 

%% GLOBAL TEST CONDITIONS /////////////////////////////////////////////////
simTime = 20;
dt = 0.05;
agentVelocities = 10;
figureSet = {'isometric','inputs'};

%% TEST CASE #1 - OFFSET ANGLE CASE ///////////////////////////////////////
CREATE THE TEST AGENTS
agentIndex{1} = agent_vectorSharing_interval();
agentIndex{2} = agent_vectorSharing_interval();
% GENERATE THE SCENARIO
objectIndex = getScenario_SDF2017offsetAngle('agents',agentIndex,'velocities',agentVelocities);
% RUN OpenMAS ON THE OBJECT SET
[DATA_1,META_1] = simulation_initialise('objects',objectIndex,'simTime',simTime,'dt',dt,'figures',figureSet);

clear agentIndex objectIndex

%% TEST CASE #2 - DIRECT COLLISION CASE ///////////////////////////////////
% CREATE THE TEST AGENTS
agentIndex{1} = agent_vectorSharing_interval();
agentIndex{2} = agent_vectorSharing_interval();
% GENERATE THE SCENARIO
objectIndex = getScenario_concentricRing('agents',agentIndex,'velocities',agentVelocities);
% RUN OpenMAS ON THE OBJECT SET
[DATA_2,META_2] = simulation_initialise('objects',objectIndex,'simTime',simTime,'dt',dt,'figures',figureSet);
                                
clear agentIndex objectIndex

%% TEST CASE #3 - THREE AGENT COCENTRIC COLLISION /////////////////////////
% CREATE THE TEST AGENTS
agentIndex{1} = agent_vectorSharing_interval();
agentIndex{2} = agent_vectorSharing_interval();
agentIndex{3} = agent_vectorSharing_interval();
% GENERATE THE SCENARIO
objectIndex = getScenario_concentricRing('agents',agentIndex,'velocities',agentVelocities);
% RUN OpenMAS ON THE OBJECT SET
[DATA_3,META_3] = simulation_initialise('objects',objectIndex,'simTime',simTime,'dt',dt,'figures',figureSet);
                                   
% RETURN TO INITIAL DIRECTORY                                 
cd(parentDir);