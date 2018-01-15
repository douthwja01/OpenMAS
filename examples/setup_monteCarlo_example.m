 %% SIMULATION SETUP (simulation_setup.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is designed to generate a multiple agent scenario from the
% agent class library inside the agents sub-directory.

% Author: James A. Douthwaite 16/05/2016

clear all; close all; 

%% ADD THE PROGRAM PATHS
addpath('environment');        
addpath('objects');  
addpath('scenarios'); 
addpath('Intlab_V7.1');   
wrkDir = pwd;

%% PREPARE PARALLEL POOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('[SETUP]\tConfiguring thread pool...\n');
% p = gcp('nocreate'); % Get current pool status
% if ~isempty(p)  
%    % A thread pool exists
%    fprintf('[SETUP]\t--> Existing thread pool found.\n'); 
%    threadObj = p;
%    %fprintf('[SETUP]\t--> closing..\n'); 
%    %delete(p) % Close any existing pool  
% else
%    fprintf('[SETUP]\t--> no active pool found.\n');
% end
% 
% % Open new pool for the monte-carlo simulations
% fprintf('[SETUP]\tOpening new thread pool...\n');
% % threadObj = parpool('IdleTimeout', 120);   % Generate thread pool, 2hour timeout
% fprintf('[SETUP]\t--> Thread pool ready.\n');

%% TEST FOR INTLAB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If intlab needs to be reloaded
try 
    test = infsup(0,1);
    clearvars test 
catch 
    IntDir = strcat(pwd,'\Intlab_V7.1\startupJD.m');
    run('startup.m'); 
    cd(wrkDir)
end 

%% 1. DEFINE THE VECTOR OF AGENT TYPES
agentNumber = 2;
for index = 1:agentNumber
    agentIndex{index} = agent_vectorSharing();
%     agentIndex{index} = agent_vectorSharing_interval();
end

%% MONTE-CARLO CYCLES SETUP 
executionDateTime = datestr(datetime('now'),' yyyy-mm-dd @ HH-MM-SS');
totalCycles = 100;
figureSet = {'fig'};
outputVector = zeros(2,totalCycles);
totalCollisions = 0; totalWaypoints = 0;
averageComputationTime = zeros(1,totalCycles);

% BEGIN SIMULATION CYCLES
for cycle = 1:totalCycles
    fprintf('[MONTE]\t MONTE-CARLO CYCLE %s \n',num2str(cycle));
    pause(rand(1));
    % GET THE INITIALISATION CONDITIONS FOR THE AGENT SET
    objectIndex = getScenario_SDF2017offsetAngle('agents',agentIndex,'noisefactor',1);
%     objectIndex = getScenario_SDF2017concentric('agents',agentIndex,'noisefactor',1);
    try
        % INITIALISE SIM WITH THE OBJECT INDEX
        [DATA,META] = simulation_initialise('objects',objectIndex,...
            'simTime',40,...
            'dt',0.25,...
            'figures',figureSet,...
        'outputPath',sprintf('MonteCarlo %s',executionDateTime));
    catch
        continue
    end
    %COLLECT OUTPUT DATA                                
    cycleCollisions = DATA.collisions;
    cycleWaypoints = DATA.waypointsAchieved;
    totalCollisions = totalCollisions + cycleCollisions;
    totalWaypoints = totalWaypoints + cycleWaypoints;
    outputVector(1,cycle) = cycleCollisions; 
    outputVector(2,cycle) = cycleWaypoints;
    % COMPUTATION TIMES
    averageComputationTime(cycle) = DATA.MEANS.global;
end
globalMean = sum(averageComputationTime)/totalCycles
totalCollisions
totalWaypoints

f = figure();
plot(linspace(1,totalCycles,totalCycles),outputVector)
grid on;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except DATA META
% load(strcat(META.outputPath,'META.mat'));
% load(strcat(META.outputPath,'EVENTS.mat'));

fprintf('..Complete.\n');