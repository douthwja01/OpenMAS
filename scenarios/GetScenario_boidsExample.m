function [ objectIndex ] = GetScenario_boidsExample(varargin)
% This function designs a simple scenario designed to demonstrate the boids
% algorithm

% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting typical concentric sphere scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',10,...
                       'agentVelocity',0,...
                       'offsetAngle',0,...
                       'waypoints',5,...
                       'waypointOrbit',10,...
                       'waypointOffsetAngle',pi/5,...
                       'waypointRadius',0.5,...
                       'noiseFactor',0,...
                       'plot',0);  
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
SBinstance = scenarioBuilder();
[scenarioConfig] = SBinstance.configurationParser(defaultConfig,varargin);

if isempty(scenarioConfig.waypointOrbit)
    scenarioConfig.waypointOrbit = scenarioConfig.agentOrbit;
end    
scenarioConfig.waypointOffsetAngle = pi + scenarioConfig.offsetAngle;        % Waypoints oppose agents
agentIndex = scenarioConfig.agents;

% DECLARE THE NUMBER OF AGENTS
agentNumber = numel(scenarioConfig.agents);
% DEFINE THE AGENT CONFIGURATIONS
agentConfig = SBinstance.regularSphere(...
    'objects',agentNumber,...
    'radius',scenarioConfig.agentOrbit,...
    'velocity',scenarioConfig.agentVelocity,...
    'zeroAngle',scenarioConfig.offsetAngle);

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
waypointConfig = SBinstance.helix(...
    'objects',scenarioConfig.waypoints,...
    'radius',scenarioConfig.waypointOrbit);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n');
% Extract the leading agent
leadAgent = agentIndex{1};

for index = 1:scenarioConfig.waypoints
    priority = scenarioConfig.waypoints - index;
    nameString = sprintf('WP%d-%s',index,leadAgent.name);
    % Construct the way-point
    waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',nameString);
    % Update the GLOBAL properties
    waypointIndex{index}.SetGLOBAL('position',waypointConfig.positions(:,index) + scenarioConfig.noiseFactor*[randn(2,1);0]);
    waypointIndex{index}.SetGLOBAL('velocity',waypointConfig.velocities(:,index));
    waypointIndex{index}.SetGLOBAL('quaternion',waypointConfig.quaternions(:,index));
    % Create the way-point association
    waypointIndex{index} = waypointIndex{index}.CreateAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end
% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    objectScenario.plotObjectIndex(objectIndex);
end
clearvars -except objectIndex % Clean-up
end