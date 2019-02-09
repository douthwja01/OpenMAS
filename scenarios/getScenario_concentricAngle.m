function [ objectIndex ] = getScenario_concentricAngle(varargin)
% This function returns the concentric scenario devised for the SDF 2017
% conference paper submission. Initialises the agents with the specific
% conditions used in second example in the paper.

fprintf('[SCENARIO]\tGetting the concentric offset-angle scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',250,...
                       'agentVelocity',18,...
                       'plot',0,...
                       'waypointOrbit',[],...
                       'waypointRadius',2,...
                       'noiseFactor',0,...
                       'angle',pi/4);
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
if isempty(scenarioConfig.waypointOrbit)
    scenarioConfig.waypointOrbit = scenarioConfig.agentOrbit;
end    
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

% GET THE SCENARIO BUILDING TOOLS
testScenario = scenarioBuilder(agentNumber);

% DEFINE THE AGENT CONFIGURATIONS
% 30mph - 13.4112m/s
% 40mph - 17.9916m/s
[ ~, agentConfig ] = testScenario.planarAngle('radius',scenarioConfig.agentOrbit,...
                                          'velocities',scenarioConfig.agentVelocity,...
                                         'offsetAngle',scenarioConfig.angle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index) + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index) + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
waypointPlanarRotation = pi;
[ ~,waypointConfig] = testScenario.planarAngle('radius',scenarioConfig.waypointOrbit,...
                                           'velocities',0,...
                                            'zeroAngle',waypointPlanarRotation,...
                                          'offsetAngle',scenarioConfig.angle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(:,index);
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
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