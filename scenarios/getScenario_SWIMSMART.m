function [ objectIndex ] = getScenario_SWIMSMART(varargin)
% This function designs the test scenario used to example the interval 
% geometric avoidance algorithm for the SWIMSMART 2017 workshop.
% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting SWIMSMART test scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',200,...
                       'velocities',18,...
                       'plot',0,...
                       'waypointRadius',5,...
                       'noiseFactor',0,...
                       'offsetAngle',pi/4);
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents
% GET THE SCENARIO BUILDING TOOLS
testScenario = scenarioBuilder(agentNumber);
% DEFINE THE AGENT CONFIGURATIONS
[ agentConfig ] = testScenario.planarRing('radius',scenarioConfig.agentOrbit,...
                                      'velocities',scenarioConfig.velocities);

%% REBUILD THE AGENT INDEX UNDER THIS CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent definitions...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);
end
% PLOT THE AGENT CONFIGURATION
% hand = agentConfig.plot();

% DEFINE THE WAYPOINT CONFIGURATIONS
waypointPlanarRotation = (18/16)*pi;
[ waypointConfig] = testScenario.planarRing('radius',scenarioConfig.agentOrbit,...
                                            'velocities',0,...
                                            'zeroAngle',waypointPlanarRotation);
                                        
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    nameString = sprintf('WP:%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(:,index);
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index});  % Create waypoint with association to agent
end
% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    testScenario.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end

% PARSE CONFIGURATION INPUTS
function [agentIndex,configFile] = parseInputs(inputSet)

% DEFAULT FILE NAME
configFile = 'scenario.mat';
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

% CHECK THERE IS SOMETHING TO SIMULATE
if isempty(agentIndex)
    error('Please provide a cell array of agent objects to simulation');
end

end