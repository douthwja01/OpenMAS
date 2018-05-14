function [ objectIndex ] = getScenario_waypoint_90Degrees(varargin)
% This scenario is designed to present a waypoint following exercise.

fprintf('[SCENARIO]\tGetting the waypoint following exercise.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'noiseFactor',0,...
                       'plot',0);  
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

%% DEFINE THE AGENT CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent definition...\n'); 
for i=1:agentNumber
    agentIndex{i}.VIRTUAL.globalPosition = [0;0;0] + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{i}.VIRTUAL.globalVelocity = [5;0;0] + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{i}.VIRTUAL.quaternion = [1;0;0;0];
end

%% DEFINE THE WAYPOINT CONFIGURATION
fprintf('[SCENARIO]\tBuilding the new scenario...\n');
waypointRadius = 0.1;
nameString = sprintf('WP:%s',agentIndex{1}.name);
waypointPositions = [ 20 20;
                       0 20;
                       0 20];

% ASSIGN AGENT GLOBAL PROPERTIES, ONE SIDE OF THE RINGS TO THE OTHER
waypointIndex{1} = waypoint('radius',waypointRadius,'name',nameString);
waypointIndex{1}.VIRTUAL.globalPosition = waypointPositions(:,1);
waypointIndex{1}.VIRTUAL.globalVelocity = zeros(3,1);
waypointIndex{1}.VIRTUAL.quaternion = [1;0;0;0];                           % Append properties from the sphereical scenario

waypointIndex{2} = waypoint('radius',waypointRadius,'name',nameString);
waypointIndex{2}.VIRTUAL.globalPosition = waypointPositions(:,2);
waypointIndex{2}.VIRTUAL.globalVelocity = zeros(3,1);
waypointIndex{2}.VIRTUAL.quaternion = [1;0;0;0];                           % Append properties from the sphereical scenario

for index = 1:numel(waypointIndex)
   waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{1},1/index);  % Create waypoint with association to agent 
end


% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
%     scenarioBuilder.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex

end