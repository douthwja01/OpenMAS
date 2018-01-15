function [ objectIndex ] = getScenario_SDF2017offsetAngle(varargin)
% This function returns the concentric scenario devised for the SDF 2017
% conference paper submission. Initialises the agents with the specific
% conditions used in second example in the paper.

fprintf('[SCENARIO]\tGetting the SDF2017 offset-angle collision scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct('file','scenario_SDF2017a.mat',...
                       'agents',[],...
                       'radius',250,...
                       'velocities',18,...
                       'plot',0,...
                       'waypointRadius',2,...
                       'noiseFactor',0,...
                       'offsetAngle',pi/4);
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

% GET THE SCENARIO BUILDING TOOLS
testScenario = scenarioBuilder(agentNumber);

% DEFINE THE AGENT CONFIGURATIONS
% 30mph - 13.4112m/s
% 40mph - 17.9916m/s
[ agentConfig ] = testScenario.planarAngle('radius',scenarioConfig.radius,...
                                           'velocities',scenarioConfig.velocities,...
                                           'offsetAngle',scenarioConfig.offsetAngle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.globalPosition = agentConfig.position(:,index) + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{index}.globalVelocity = agentConfig.velocity(:,index) + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{index}.quaternion = agentConfig.quaternion(:,index);
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
waypointPlanarRotation = pi;
[ waypointConfig] = testScenario.planarAngle('radius',scenarioConfig.radius,'velocities',0,...
                                             'zeroAngle',waypointPlanarRotation,...
                                             'offsetAngle',scenarioConfig.offsetAngle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    nameString = sprintf('WP:%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end
% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    scenarioBuilder.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end

% PARSE CONFIGURATION INPUTS
function [agentIndex,noiseFactor,configFile] = parseInputs(inputSet)

% DEFAULT FILE NAME
configFile = 'scenario.mat';
noiseFactor = 0;
agentIndex = {};

% VISABILTY REDUCTION RATIO (IN FUTURE PROPORTIONAL TO ISO PROPERTIES)
tmp = strncmpi(inputSet,'filename',4);
if any(tmp)
    configFile = inputSet{find(tmp) + 1};
end

tmp = strncmpi(inputSet,'agents',5)|strncmpi(inputSet,'objects',6); 
if any(tmp)
    agentIndex = inputSet{find(tmp) + 1}; 
    assert(iscell(agentIndex) == 1,'Figure set must be specified as a cell array of string names.');
end

tmp = strncmpi(inputSet,'noisefactor',5); % |strncmpi(inputSet,'objects',6); 
if any(tmp)
    noiseFactor = inputSet{find(tmp) + 1}; 
    assert(isnumeric(noiseFactor) == 1,'Noise factor must be a numeric factor.');
end

% CHECK THERE IS SOMETHING TO SIMULATE
if isempty(agentIndex)
    error('Please provide a cell array of agent objects to simulation');
end

end