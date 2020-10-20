%% OpenMAS THESIS TEST (setup_thesis.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the comparison of the VO methods for the exampled
% methods in the thesis document.

% Author: James A. Douthwaite 26/03/18

% ADD THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('toolboxes');
addpath('scenarios'); 

% CODE DEBUG COMMANDS %
% profile on
% profile viewer 

%If intlab needs to be reloaded
%IntLab();   % Load Intlab

% ///////////////////////// SCENARIO PARAMETERS ///////////////////////////
fprintf('[SETUP]\tInitialising the comparative example script.\n');
SCENARIO = struct();
SCENARIO.agentOrbit    = 10;       
SCENARIO.agentVelocity = 0;      
SCENARIO.waypointOrbit = 15;               % Parameters used to define the scenario
SCENARIO.offsetAngle = pi/2;
SCENARIO.noiseSigma = 0.1;
SCENARIO.plot = false;

% ///////////////////////// AGENT CONFIGURATION ///////////////////////////
SCENARIO.case = 1;                         % The scenario number (1-3)
SCENARIO.agentNumber = 2:5;                % Total agent number
SCENARIO.algorithms  = {...
    'agent_2D_VO',...
    'agent_2D_RVO',...
    'agent_2D_HRVO',...
    'agent_2D_ORCA'};                      % The selected agents within the '/objects' folder
% SCENARIO.algorithms = {...
%     'agent_2D_IA',...
%     'agent_IA'};                                                            

% //////////////////////// OMAS CONFIGURATION /////////////////////////////
OMAS = struct();
OMAS.maxDuration = 60;    
OMAS.timeStep = 0.25;      
OMAS.verbosity = 0;                        % Parameters defining the OMAS configuration within a give cycle. 

% SIMULATION PARAMETERS
% sim_figureSet ={'all'};
OMAS.figureSet = {'none'};
% sim_figureSet = {'plan','closest','inputs','plan'}; 

% ///////////////////// MONTE-CARLO CONFIGURATION /////////////////////////
[~, userdir]  = system('echo %USERPROFILE%'); % Get desktop path
MC = struct();
MC.outputPath = strcat(userdir,'\desktop\openmas-data'); 
MC.cycles     = 10;
MC.threadPool = false;

% ////////////////////// PROGRESSIVE SIMULATIONS //////////////////////////
MC.index = {};
for i = 1:numel(SCENARIO.agentNumber)
    % Move through the experiments with respect to each algorithm.
    for j = 1:numel(SCENARIO.algorithms)
        % //////// DEFINE THE AGENT CONFIGURATION /////////       
        fprintf('[SETUP]\tDefining the agent configuration.\n');
        
        % Create the agent set
        agentIndex = cell(1,SCENARIO.agentNumber(i));                             % Agent container
        for k = 1:SCENARIO.agentNumber(i) 
            agentIndex{1,k} = eval(SCENARIO.algorithms{j});                    % Assign agents
        end
        
        % //////// GET THE SCENARIO CONFIGURATION /////////
        objectIndex = DefineScenario(SCENARIO,agentIndex);
        % /////////////////////////////////////////////////

        % Assemble array
        [MC.index] = vertcatPadded(MC.index,objectIndex);
    end
end
% Create a Monte-Carlo analysis instance
MCinstance = OMAS_monteCarlo(...
    'objects',MC.index,...
    'cycles',MC.cycles,...
    'isParallel',MC.threadPool,...
    'duration',OMAS.maxDuration,...
    'dt',OMAS.timeStep,...
    'figures',OMAS.figureSet,...
    'directory',MC.outputPath,...
    'positionSigma',0.2,...
    'velocitySigma',0.2); 

% Begin evaluating the cycles
MCinstance = MCinstance.EvaluateAllCycles();

% /////////////////// SCENARIO CONFIGURATION //////////////////////////
function [objectIndex] = DefineScenario(SCENARIO,agentIndex)
switch SCENARIO.case
    case 1
        objectIndex = GetScenario_concentricRing(...                   % Initialise the objects in the UKACC scenario
            'agents',agentIndex,...
            'agentOrbit',SCENARIO.agentOrbit,...
            'agentVelocity',SCENARIO.agentVelocity,...
            'waypointOrbit',SCENARIO.waypointOrbit,...
            'noiseFactor',SCENARIO.noiseSigma,...
            'plot',SCENARIO.plot);
    case 2
        objectIndex = GetScenario_concentricAngle(...
            'agents',agentIndex,...
            'agentOrbit',SCENARIO.agentOrbit,...
            'agentVelocity',SCENARIO.agentVelocity,...
            'waypointOrbit',SCENARIO.waypointOrbit,...
            'angle',SCENARIO.offsetAngle,...
            'plot',SCENARIO.plot);
    case 3
        objectIndex = GetScenario_twoLines(...
            'agents',agentIndex,...
            'agentSeparation',SCENARIO.agentOrbit,...
            'agentVelocity',SCENARIO.agentVelocity,...
            'padding',5,...
            'waypointSeparation',SCENARIO.waypointOrbit,...
            'noiseFactor',SCENARIO.noiseSigma,...
            'plot',SCENARIO.plot);
    otherwise
        error('Scenario number not recognised');
end
% Reformat the array
objectIndex = reshape(objectIndex,[1,numel(objectIndex)]);
end