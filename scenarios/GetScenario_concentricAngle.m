function [ objectIndex ] = GetScenario_concentricAngle(varargin)
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
                       'waypointRadius',0.5,...
                       'noiseFactor',0,...
                       'angle',pi/4);
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% Check
if isempty(inputConfig.waypointOrbit)
    inputConfig.waypointOrbit = inputConfig.agentOrbit;
end    
agentIndex = inputConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

% DEFINE THE AGENT CONFIGURATIONS
% 30mph - 13.4112m/s
% 40mph - 17.9916m/s
agentConfig = SBinstance.planarAngle(...
    'objects',agentNumber,...
    'radius',inputConfig.agentOrbit,...
    'velocity',inputConfig.agentVelocity,...
    'offsetAngle',inputConfig.angle);

%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocities(:,index) + inputConfig.noiseFactor*randn(3,1);
    agentIndex{index}.VIRTUAL.quaternion     = agentConfig.quaternions(:,index);
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
waypointPlanarRotation = pi;
waypointConfig = SBinstance.planarAngle(...
    'objects',agentNumber,...
    'radius',inputConfig.waypointOrbit,...
    'velocities',0,...
    'zeroAngle',waypointPlanarRotation,...
    'offsetAngle',inputConfig.angle);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:agentNumber
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.positions(:,index);
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