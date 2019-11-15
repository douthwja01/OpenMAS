function [ objectIndex ] = GetScenario_waypointCurve(varargin)
% This scenario is designed to present a waypoint following exercise.

fprintf('[SCENARIO]\tGetting the waypoint following exercise.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentVelocity',0,...
                       'noiseFactor',0,...
                       'waypoints',3,...
                       'waypointRadius',0.1,...
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
    agentIndex{i}.SetGLOBAL('position',[0;0;0] + inputConfig.noiseFactor*randn(3,1));
    agentIndex{i}.SetGLOBAL('velocity',[inputConfig.agentVelocity;0;0] + inputConfig.noiseFactor*randn(3,1));
    agentIndex{i}.SetGLOBAL('quaternion',[1;0;0;0]);
end

%% DEFINE THE WAYPOINT CONFIGURATION
fprintf('[SCENARIO]\tBuilding the new scenario...\n');
waypointDefiningRadius = 20;
angleOffset = -pi/2;                                % Align the first waypoint to be directly infront of the agent
waypointConfig = SBinstance.planarAngle(...
    'objects',agentNumber,...
    'radius',waypointDefiningRadius,...
    'pointA',[waypointDefiningRadius;-waypointDefiningRadius;-1],...
    'pointB',[waypointDefiningRadius;-waypointDefiningRadius;0],...
    'velocities',0,...
    'zeroAngle',angleOffset);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:inputConfig.waypoints
    nameString = sprintf('WP-%s',agentIndex{1}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.SetGLOBAL('position',waypointConfig.positions(:,index));
    waypointIndex{index}.SetGLOBAL('velocity',waypointConfig.velocities(:,index));
    waypointIndex{index}.SetGLOBAL('quaternion',waypointConfig.quaternions(:,index));
    % Assign waypoint to agent with priority
    waypointIndex{index} = waypointIndex{index}.CreateAgentAssociation(agentIndex{1},1/index);  % Create waypoint with association to agent
end

% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
clearvars -except objectIndex % Clean-up
end