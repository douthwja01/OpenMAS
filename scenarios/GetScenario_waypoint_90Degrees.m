function [ objectIndex ] = GetScenario_waypoint_90Degrees(varargin)
% This scenario is designed to present a waypoint following exercise.

fprintf('[SCENARIO]\tGetting the waypoint following exercise.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentVelocity',0,...
                       'noiseFactor',0,...
                       'plot',0);  

% Instanciate the scenario builder
SBinstance = scenarioBuilder();                  
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex = inputConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

%% DEFINE THE AGENT CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent definition...\n'); 
for i=1:agentNumber
    agentIndex{i}.SetGLOBAL('globalPosition',[0;0;0] + inputConfig.noiseFactor*randn(3,1));
    agentIndex{i}.SetGLOBAL('globalVelocity',[inputConfig.agentVelocity;0;0] + inputConfig.noiseFactor*randn(3,1));
    agentIndex{i}.SetGLOBAL('quaternion',[1;0;0;0]);
end

%% DEFINE THE WAYPOINT CONFIGURATION
fprintf('[SCENARIO]\tBuilding the new scenario...\n');
waypointRadius = 0.1;
nameString = sprintf('WP-%s',agentIndex{1}.name);
waypointPositions = [ 20 20;
                       0 20;
                       0 20];

% ASSIGN AGENT GLOBAL PROPERTIES, ONE SIDE OF THE RINGS TO THE OTHER
waypointIndex{1} = waypoint('radius',waypointRadius,'name',nameString);
waypointIndex{1}.SetGLOBAL('globalPosition',waypointPositions(:,1));
waypointIndex{1}.SetGLOBAL('globalVelocity',zeros(3,1));
waypointIndex{1}.SetGLOBAL('quaternion',[1;0;0;0]);                        % Append properties from the sphereical scenario

waypointIndex{2} = waypoint('radius',waypointRadius,'name',nameString);
waypointIndex{2}.SetGLOBAL('globalPosition',waypointPositions(:,2));
waypointIndex{2}.SetGLOBAL('globalVelocity',zeros(3,1));
waypointIndex{2}.SetGLOBAL('quaternion',[1;0;0;0]);                        % Append properties from the sphereical scenario

for index = 1:numel(waypointIndex)
   waypointIndex{index} = waypointIndex{index}.CreateAgentAssociation(agentIndex{1},1/index);  % Create waypoint with association to agent 
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