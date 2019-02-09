function [ objectIndex ] = getScenario_obstacleTrack(varargin)
% This function is designed to generate a simple cube of obstacles 

fprintf('[SCENARIO]\tGetting Simple thre obstacle, one waypoint scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...   
                       'agentVelocity',0,...
                       'obstacles',[],...
                       'obstacleRadius',1,...
                       'waypoints',[],...
                       'waypointRadius',0.5,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[config] = scenarioBuilder.configurationParser(defaultConfig,varargin);

% INPUT SANITY CHECK
assert(isnumeric(config.agentVelocity) && numel(config.obstacles) == 1,'Agent velocity must be a scalar.');
assert(isnumeric(config.obstacleRadius) && numel(config.obstacleRadius) == 1,'Obstacle radius must be a scalar.');
assert(isnumeric(config.waypointRadius) && numel(config.waypointRadius) == 1,'Waypoint radius must be a scalar.');
assert(iscell(config.agents) && numel(config.agents) == 3,'This scenario is intended to use three obstacles only.');
assert(iscell(config.obstacles) && numel(config.obstacles) == 3,'This scenario is intended to use three obstacles only.');

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentIndex = config.agents;  
% Place the agent in the global axis 
agentIndex{1}.VIRTUAL.radius = config.agentRadius;
% APPLY GLOBAL STATE VARIABLES
agentIndex{1}.VIRTUAL.globalPosition = [0;0;0];
agentIndex{1}.VIRTUAL.globalVelocity = [1;0;0]*config.agentVelocity;
agentIndex{1}.VIRTUAL.quaternion = [1;0;0;0];

%% ////////////////// BUILD THE OBSTACLE GLOBAL STATES ////////////////////
obstacleNumber = numel(config.obstacles); 
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
% POSITION COEFFICIENTS
diffSpacing = 2*config.obstacleRadius;
diagSpacing = 4*config.obstacleRadius;
xyVector = [0.5;-0.5];

for index = 1:obstacleNumber
    obstacleIndex{index} = config.obstacles{index};                   % Get the agents from the input structure
    obstacleIndex{index}.name = sprintf('OB-%s',config.obstacles{index}.name);
    % APPLY GENERIC GLOBAL STATE VARIABLES
    obstacleIndex{index}.VIRTUAL.globalVelocity = [0;0;0];
    obstacleIndex{index}.VIRTUAL.quaternion = [1;0;0;0];  % Append properties from the sphereical scenario
    
    % GENERATE POSITIONS
%     position = ((diagSpacing*index)+diffSpacing*xyVector);
    obstacleIndex{index}.VIRTUAL.globalPosition = [((diagSpacing*index)+diffSpacing*xyVector);0];
    xyVector = -xyVector; % OSCILLATE THE Y POSITIONS
end


%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
waypointIndex = cell(1);

nameString = sprintf('WP-%s',agentIndex{1}.name);
waypointIndex{1} = waypoint('radius',config.waypointRadius,'priority',1,'name',nameString);
% APPLY GLOBAL STATE VARIABLES
waypointIndex{1}.VIRTUAL.globalPosition = [(diagSpacing*(index+1));(diagSpacing*(index+1));0];
waypointIndex{1}.VIRTUAL.globalVelocity = [0;0;0];
waypointIndex{1}.VIRTUAL.quaternion = [1;0;0;0];
waypointIndex{1} = waypointIndex{1}.createAgentAssociation(agentIndex{1});  % Create waypoint with association to agent

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex;obstacleIndex;waypointIndex]; 
% SAVE THE FILE
save(config.file,'objectIndex');
% PLOT THE SCENE
if config.plot
    SBInstance = scenarioBuilder(numel(objectIndex));
    SBInstance.plotObjectIndex(objectIndex);                          % Plot the object index
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end