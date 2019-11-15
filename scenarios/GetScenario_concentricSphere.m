function [ objectIndex ] = GetScenario_concentricSphere(varargin)
% This function designs a typical three agent, three waypoint collision
% scenario.

% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting typical concentric sphere scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentOrbit',10,...
    'agentVelocity',0,...
    'offsetAngle',0,...
    'waypointOrbit',[],...
    'waypointOffsetAngle',[],...
    'waypointRadius',0.5,...
    'noiseFactor',0,...
    'plot',false);  

% Instantiae the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% Check
if isempty(inputConfig.waypointOrbit)
    inputConfig.waypointOrbit = inputConfig.agentOrbit;
end    
inputConfig.waypointOffsetAngle = pi + inputConfig.offsetAngle;        % Waypoints oppose agents
agentIndex = inputConfig.agents;

% DECLARE THE NUMBER OF AGENTS
agentNumber = numel(inputConfig.agents);

% DEFINE THE AGENT CONFIGURATIONS
agentConfig = SBinstance.regularSphere(...
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
% Define the object-index
objectIndex = horzcat(agentIndex); %,waypointIndex);
% Plot the scene
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
clearvars -except objectIndex % Clean-up
end