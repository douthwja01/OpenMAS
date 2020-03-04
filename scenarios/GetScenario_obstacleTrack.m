function [ objectIndex ] = GetScenario_obstacleTrack(varargin)
% This function is designed to generate a simple cube of obstacles 

fprintf('[SCENARIO]\tGetting Simple thre obstacle, one waypoint scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...   
    'agentVelocity',0,...
    'obstacles',[],...
    'obstacleRadius',1,...
    'waypoints',[],...
    'waypointRadius',0.5,...
    'plot',0);
                   
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);

% INPUT SANITY CHECK
assert(isnumeric(inputConfig.agentVelocity) && numel(inputConfig.agentVelocity) == 1,'Agent velocity must be a scalar.');
assert(isnumeric(inputConfig.obstacleRadius) && numel(inputConfig.obstacleRadius) == 1,'Obstacle radius must be a scalar.');
assert(isnumeric(inputConfig.waypointRadius) && numel(inputConfig.waypointRadius) == 1,'Waypoint radius must be a scalar.');
% Check agents
assert(iscell(inputConfig.agents) && numel(inputConfig.agents) == 1,'This scenario is intended to use one agent only.');
% Check obstacles
assert(iscell(inputConfig.obstacles) && numel(inputConfig.obstacles) == 3,'This scenario is intended to use three obstacles only.');
assert(iscell(inputConfig.obstacles),'This scenario is intended to use three obstacles only.');

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentIndex = inputConfig.agents;  
% APPLY GLOBAL STATE VARIABLES
agentIndex{1}.SetGLOBAL('position',[0;0;0]);
agentIndex{1}.SetGLOBAL('velocity',[1;0;0]*inputConfig.agentVelocity);

%% ////////////////// BUILD THE OBSTACLE GLOBAL STATES ////////////////////
obstacleNumber = numel(inputConfig.obstacles); 
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
% POSITION COEFFICIENTS
diffSpacing = 2*inputConfig.obstacleRadius;
diagSpacing = 4*inputConfig.obstacleRadius;
xyVector = [0.5;-0.5];

for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};                   % Get the agents from the input structure
    obstacleIndex{index}.name = sprintf('OB-%s',inputConfig.obstacles{index}.name);
    % APPLY GENERIC GLOBAL STATE VARIABLES
    obstacleIndex{index}.SetGLOBAL('position',[((diagSpacing*index)+diffSpacing*xyVector);0]);
    obstacleIndex{index}.SetGLOBAL('velocity',[0;0;0]);
    obstacleIndex{index}.SetGLOBAL('quaternion',[1;0;0;0]);  % Append properties from the sphereical scenario
    xyVector = -xyVector; % Define the opposing vector
end

%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
waypointIndex = cell(1);

nameString = sprintf('WP-%s',agentIndex{1}.name);
waypointIndex{1} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
% APPLY GLOBAL STATE VARIABLES
waypointIndex{1}.SetGLOBAL('position',[(diagSpacing*(index+1));(diagSpacing*(index+1));0];
waypointIndex{1}.SetGLOBAL('velocity',[0;0;0]);
waypointIndex{1}.SetGLOBAL('quaternion',[1;0;0;0]);
waypointIndex{1} = waypointIndex{1}.CreateAgentAssociation(agentIndex{1});  % Create waypoint with association to agent

%% UPDATE THE AGENTS INITIAL HEADING
% Design the agents initial global heading
headingVector = waypointIndex{1}.GetGLOBAL('position') - agentIndex{1}.GetGLOBAL('position');
headingVector = headingVector/norm(headingVector);
q = OMAS_geometry.vectorsToQuaternion(headingVector,[1;0;0]);
agentIndex{1}.GetGLOBAL('quaternion',q);

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex,obstacleIndex,waypointIndex]; 
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);                          % Plot the object index
end
clearvars -except objectIndex % Clean-up 
end