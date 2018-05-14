%% THREE OBSTACLE- THREE AGENT SCENARIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a three obstacle, three agent scenario. This
% scenario is designed to evaluate a SWARMS ability to converge on a
% desired target formation in the prescence of a obstacle avoidance
% algorithm.

% The scenario consists of the following:
% - 3x generic agents
% - 3x waypoints

function [ objectIndex ] = getScenario_threeObstacles_shiyu(varargin)

fprintf('[SCENARIO]\tDESIGNING FORMATION-AVOIDANCE THREE OBSTACLE TEST CASE.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct('obstacleSize',3,...
                       'obstacleOrbit',10,...
                       'agentOrbit',20,...
                       'velocities',0.5,...
                       'waypointSize',0.1,...
                       'waypointOrbit',2,...
                       'obstacles',[],...
                       'agents',[],...
                       'plot',0);
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;

for i=1:3
    obstacleIndex{i} = obstacle('size',scenarioConfig.obstacleSize);               % Create the obstacle set
end
% COUNT OBSTACLES
obstacleNumber = numel(obstacleIndex);
% GET THE SCENARIO BUILDER INSTANCE 
testScenario = scenarioBuilder(obstacleNumber);
% DEFINE THE AGENT CONFIGURATIONS
[ obstacleConfig ] = testScenario.planarRing('radius',scenarioConfig.obstacleOrbit,'velocities',0);
% MOVE THROUGH THE OBSTACLES AND INITIALISE WITH GLOBAL PROPERTIES
for index = 1:obstacleNumber
    % APPLY GLOBAL STATE VARIABLES
    obstacleIndex{index}.VIRTUAL.globalPosition = obstacleConfig.position(index,:)';% + randn(3,1); % Apply some noise to the obstacle global positions
    obstacleIndex{index}.VIRTUAL.globalVelocity = obstacleConfig.velocity(index,:)';
    obstacleIndex{index}.VIRTUAL.quaternion = obstacleConfig.quaternion(index,:)';
end

%% DESIGN AGENT CONFIGURATION
agentNumber = numel(agentIndex);
% GET THE SCENARIO BUILDER INSTANCE 
testScenario = scenarioBuilder(agentNumber);
% DEFINE THE AGENT CONFIGURATIONS
[ agentConfig ] = testScenario.planarRing('radius',scenarioConfig.agentOrbit,'velocities',scenarioConfig.velocities);
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(index,:)' + randn(3,1); % Apply some noise to the agent global positions
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(index,:)';
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(index,:)';
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
angleOffset = pi;
[ waypointConfig] = testScenario.planarRing('radius',scenarioConfig.waypointOrbit,'velocities',0,...
                                            'zeroAngle',angleOffset);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
for index = 1:agentNumber
    nameString = sprintf('WP:%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointSize,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(index,:)';
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(index,:)';
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(index,:)';
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end

% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,obstacleIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    scenarioBuilder.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end