function [ objectIndex ] = getScenario_obstacleTrack(varargin)


fprintf('[SCENARIO]\tGetting Simple thre obstacle, one waypoint scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...   
                       'agentVelocity',0,...
                       'agentRadius',0.5,...    % Width of 1m
                       'obstacles',3,...
                       'obstacleRadius',1,...
                       'waypoints',[],...
                       'waypointRadius',0.5,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex = inputConfig.agents;                                           % Declare the agent set
agentNumber = numel(agentIndex);                                           % Declare the number of agents
assert(agentNumber == 1,'This scenario is intended to test an individual agent.');
% GENERATE OBSTACLE OBJECTS
if isnumeric(inputConfig.obstacles)
    obstacleIndex = cell(inputConfig.obstacles,1);
    for index = 1:inputConfig.obstacles
       obstacleIndex{index} = obstacle();
    end
    inputConfig.obstacles = obstacleIndex;
end
% DECLARE THE NUMBER OF OBSTACLES
obstacleNumber = numel(inputConfig.obstacles);  
assert(obstacleNumber == 3,'This scenario is intended to use three obstacles only.');

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
% Place the agent in the global axis 
agentIndex{1}.VIRTUAL.radius = inputConfig.agentRadius;
% APPLY GLOBAL STATE VARIABLES
agentIndex{1}.VIRTUAL.globalPosition = [0;0;0];
agentIndex{1}.VIRTUAL.globalVelocity = [1;0;0]*inputConfig.agentVelocity;
agentIndex{1}.VIRTUAL.quaternion = [1;0;0;0];

%% ////////////////// BUILD THE OBSTACLE GLOBAL STATES ////////////////////
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
% POSITION COEFFICIENTS
diffSpacing = 2*inputConfig.obstacleRadius;
diagSpacing = 4*inputConfig.obstacleRadius;
xyVector = [0.5;-0.5];

for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};                   % Get the agents from the input structure
    obstacleIndex{index}.name = sprintf('OB-%s',inputConfig.obstacles{index}.name);
    obstacleIndex{index}.VIRTUAL.radius = inputConfig.obstacleRadius;
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
waypointIndex{1} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
% APPLY GLOBAL STATE VARIABLES
waypointIndex{1}.VIRTUAL.globalPosition = [(diagSpacing*(index+1));(diagSpacing*(index+1));0];
waypointIndex{1}.VIRTUAL.globalVelocity = [0;0;0];
waypointIndex{1}.VIRTUAL.quaternion = [1;0;0;0];
waypointIndex{1} = waypointIndex{1}.createAgentAssociation(agentIndex{1});  % Create waypoint with association to agent

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex;obstacleIndex;waypointIndex]; 
% SAVE THE FILE
save(inputConfig.file,'objectIndex');
% PLOT THE SCENE
if inputConfig.plot
    SBInstance = scenarioBuilder(numel(objectIndex));
    SBInstance.plotObjectIndex(objectIndex);                          % Plot the object index
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end