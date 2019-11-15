function [ objectIndex ] = GetScenario_concentricRing(varargin)
% This function designs a typical three agent, three waypoint collision
% scenario.

% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting typical concentric agent scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentOrbit',10,...
    'agentVelocity',0,...
    'offsetAngle',0,...
    'waypointOrbit',[],...
    'waypointOffsetAngle',[],...
    'waypointRadius',0.1,...
    'noiseFactor',0,...
    'plot',false);  

% Instanciate the scenario builder
SBinstance = scenarioBuilder();               
% Parse user inputs
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% Check
if isempty(inputConfig.waypointOrbit)
    inputConfig.waypointOrbit = inputConfig.agentOrbit;
end    
inputConfig.waypointOffsetAngle = pi + inputConfig.offsetAngle;        % Waypoints oppose agents
agentIndex = inputConfig.agents;
% Declare the agent
agentNumber = numel(inputConfig.agents);

% Define the agent configuration
agentConfig = SBinstance.planarRing(...
    'objects',agentNumber,...
    'radius',inputConfig.agentOrbit,...
    'velocity',inputConfig.agentVelocity,...
    'zeroAngle',inputConfig.offsetAngle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % Update the GLOBAL properties
    agentIndex{index}.SetGLOBAL('position',agentConfig.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]); % 2D PERTURBATION
    agentIndex{index}.SetGLOBAL('velocity',agentConfig.velocities(:,index));
    agentIndex{index}.SetGLOBAL('quaternion',agentConfig.quaternions(:,index));
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
waypointConfig = SBinstance.planarRing(...
    'objects',agentNumber,...
    'radius',inputConfig.waypointOrbit,...
    'velocities',0,...
    'zeroAngle',inputConfig.waypointOffsetAngle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    % Create a way-point
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'name',sprintf('WP-%s',agentIndex{index}.name));
    % Update the GLOBAL properties
    waypointIndex{index}.SetGLOBAL('position',waypointConfig.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]);
    waypointIndex{index}.SetGLOBAL('velocity',waypointConfig.velocities(:,index));
    waypointIndex{index}.SetGLOBAL('quaternion',waypointConfig.quaternions(:,index));
    % Create the way-point association
    waypointIndex{index}.CreateAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end
% Extend the object-index
objectIndex = horzcat(agentIndex,waypointIndex);
% Plot the scene
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% Clean-up
clearvars -except objectIndex
end