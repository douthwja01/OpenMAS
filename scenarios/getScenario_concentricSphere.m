function [ objectIndex ] = getScenario_concentricSphere(varargin)
% This function designs a typical three agent, three waypoint collision
% scenario.

% The scenario consists of the following:
% - 3x agent_vectorSharing_interval agents.
% - 3x waypoints
% The agents are positioned in a ring, with the waypoints at the apposing
% side.

fprintf('[SCENARIO]\tGetting typical cocentric sphere scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',10,...
                       'velocities',0,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
                   
% GET SCENARIO BUILDER INSTANCE 
objectNumber = numel(scenarioConfig.agents);
testScenario = scenarioBuilder(objectNumber);
errorFlag = 1;
while errorFlag ~= 0        %% <-- Temporary error fix
    try  
        [agentConfig] = testScenario.randomSphere('radius',scenarioConfig.agentOrbit,'velocity',scenarioConfig.velocities);
        errorFlag = 0;
    catch
    end
end
    
%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
agentIndex = cell(objectNumber,1);
for index = 1:objectNumber
    agentIndex{index} = scenarioConfig.agents{index};
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);
end

%% DEFINE WAYPOINTS AND ASSIGN GLOBAL PARAMETERS
% DEFINE THE FIRST WAYPOINT SET
% [waypointConfig] = testScenario.planarRing('radius',1.5*scenarioConfig.radius,'velocities',0,...
%                                             'zeroAngle',scenarioConfig.waypointOffsetAngle);
% 
% % MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
% fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
% for index = 1:agentNumber
%     nameString = sprintf('WP:%s',agentIndex{index}.name);
%     waypointIndex{index} = waypoint('radius',scenarioConfig.waypointRadius,'name',nameString);
%     % APPLY GLOBAL STATE VARIABLES
%     waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(index,:)';
%     waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(index,:)';
%     waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(index,:)';
%     waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
% end
% 
% % BUILD THE COLLECTIVE OBJECT INDEX
% objectIndex = horzcat(agentIndex,waypointIndex);

objectIndex = agentIndex;

% % PLOT THE SCENE
if scenarioConfig.plot
    testScenario.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end