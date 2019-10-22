function [ objectIndex ] = GetScenario_concentricRing(varargin)
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
                       'agentOrbit',10,...
                       'agentVelocity',0,...
                       'offsetAngle',0,...
                       'waypointOrbit',[],...
                       'waypointOffsetAngle',[],...
                       'waypointRadius',0.1,...
                       'noiseFactor',0,...
                       'plot',0);  

% GET THE SCENARIO BUILDING TOOLS
SBinstance = scenarioBuilder();               
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);

if isempty(inputConfig.waypointOrbit)
    inputConfig.waypointOrbit = inputConfig.agentOrbit;
end    
inputConfig.waypointOffsetAngle = pi + inputConfig.offsetAngle;        % Waypoints oppose agents
agentIndex = inputConfig.agents;

% DECLARE THE NUMBER OF AGENTS
agentNumber = numel(inputConfig.agents);

% DEFINE THE AGENT CONFIGURATIONS
agentConfig = SBinstance.planarRing(...
    'objects',agentNumber,...
    'radius',inputConfig.agentOrbit,...
    'velocity',inputConfig.agentVelocity,...
    'zeroAngle',inputConfig.offsetAngle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]; % 2D PERTURBATION
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocities(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternions(:,index);
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
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1);
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocities(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternions(:,index);
    waypointIndex{index} = waypointIndex{index}.CreateAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end
% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end