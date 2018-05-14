function [ objectIndex ] = getScenario_waypointCurve(varargin)
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
waypointNumber = 3;
testScenario = scenarioBuilder(waypointNumber);
angleOffset = -(10/9)*pi;
[ waypointConfig] = testScenario.planarAngle('radius',20,...
                                             'pointA',[20;-20;-1],...
                                             'pointB',[20;-20;0],...
                                             'velocities',0,...
                                             'zeroAngle',angleOffset);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:waypointNumber
    waypointIndex{index} = waypoint('radius',0.1);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(:,index);
    waypointPriority = 1/index;
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{1},waypointPriority);  % Create waypoint with association to agent
end

% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    scenarioBuilder.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex

end