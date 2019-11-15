function [ objectIndex ] = GetScenario_concentricAngle(varargin)
% This function returns the concentric scenario devised for the SDF 2017
% conference paper submission. Initialises the agents with the specific
% conditions used in second example in the paper.

fprintf('[SCENARIO]\tGetting the concentric offset-angle scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentOrbit',250,...
    'agentVelocity',18,...
    'waypointOrbit',[],...
    'waypointRadius',0.5,...
    'noiseFactor',0,...
    'angle',pi/4,...
    'plot',false);

% Instantiate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[scenarioConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% Check
if isempty(scenarioConfig.waypointOrbit)
    scenarioConfig.waypointOrbit = scenarioConfig.agentOrbit;
end    
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

% Define the agent configuration
% 30mph - 13.4112m/s
% 40mph - 17.9916m/s
agentConfig = SBinstance.planarAngle(...
    'objects',agentNumber,...
    'radius',scenarioConfig.agentOrbit,...
    'velocity',scenarioConfig.agentVelocity,...
    'offsetAngle',scenarioConfig.angle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % Update the GLOBAL properties
    agentIndex{index}.SetGLOBAL('position',agentConfig.positions(:,index) + scenarioConfig.noiseFactor*[randn(2,1);0]); % 2D PERTURBATION
    agentIndex{index}.SetGLOBAL('velocity',agentConfig.velocities(:,index));
    agentIndex{index}.SetGLOBAL('quaternion',agentConfig.quaternions(:,index));
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
waypointAgentOffset = pi;
waypointConfig = SBinstance.planarAngle(...
    'objects',agentNumber,...
    'radius',scenarioConfig.waypointOrbit,...
    'velocities',0,...
    'zeroAngle',waypointAgentOffset,...
    'offsetAngle',scenarioConfig.angle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',sprintf('WP-%s',agentIndex{index}.name));
    % Update the GLOBAL properties
    waypointIndex{index}.SetGLOBAL('position',waypointConfig.positions(:,index) + scenarioConfig.noiseFactor*[randn(2,1);0]);
    waypointIndex{index}.SetGLOBAL('velocity',waypointConfig.velocities(:,index));
    waypointIndex{index}.SetGLOBAL('quaternion',waypointConfig.quaternions(:,index));
    % Create the way-point association
    waypointIndex{index} = waypointIndex{index}.CreateAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end
% Extend the object-index
objectIndex = horzcat(agentIndex,waypointIndex);
% Plot the scene
if scenarioConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% Clean-up
clearvars -except objectIndex
end