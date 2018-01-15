function [ objectIndex ] = getScenario_concentricRing(varargin)
% This function designs a typical three agent, three waypoint collision
% scenario.

% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting typical Cocentric test scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'radius',10,...
                       'waypointOffsetAngle',pi,...
                       'waypointRadius',2,...
                       'velocities',0,...
                       'plot',0);  
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;

% DECLARE THE NUMBER OF AGENTS
agentNumber = numel(scenarioConfig.agents);
% GET THE SCENARIO BUILDING TOOLS
testScenario = scenarioBuilder(agentNumber);
% DEFINE THE AGENT CONFIGURATIONS
[ agentConfig ] = testScenario.planarRing('radius',scenarioConfig.radius,...
                                          'velocities',scenarioConfig.velocities);
  
%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.quaternion = agentConfig.quaternion(:,index);
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
[waypointConfig] = testScenario.planarRing('radius',1.5*scenarioConfig.radius,'velocities',0,...
                                            'zeroAngle',scenarioConfig.waypointOffsetAngle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    nameString = sprintf('WP:%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.quaternion = waypointConfig.quaternion(:,index);
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