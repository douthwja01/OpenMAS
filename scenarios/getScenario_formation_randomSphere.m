function [ objectIndex ] = getScenario_formation_randomSphere(varargin)
% This function designs a scenario where the agents are randomly
% posistioned around a outer sphere, where the obstacles are positioned
% regularly around an inner sphere.

fprintf('[SCENARIO]\tGetting typical cocentric sphere scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',20,...
                       'obstacles',6,...                                   % Minimum required to demonstrate 3D
                       'obstacleRadius',2,...
                       'obstacleOrbit',10,...
                       'velocities',0,...
                       'adjacencyMatrix',[],...                            % The globally specified adjacency matrix
                       'plot',0);
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
                   
% THIS AGENT REQUIRES ONLY TWO AGENTS, CONFIRM
agentNumber = numel(inputConfig.agents);                                   % The number of agents

if isnumeric(inputConfig.obstacles) || isempty(inputConfig.obstacles)
    obstacleSet = cell(inputConfig.obstacles,1);
    for index = 1:inputConfig.obstacles
       obstacleSet{index} = obstacle();
    end
    inputConfig.obstacles = obstacleSet;
end  
obstacleNumber = numel(inputConfig.obstacles);                             % The number of obstacles

% DESIGN THE DESIRED SEPERATION MATRIX (ADJACENCY MATRIX)
% The adjacency matrix is indexed by objectID in scenarioConfig.adjacencyMatrix
% OTHERWISE ASSIGN DEFAULT ADJACENCY MATRIX
if isempty(inputConfig.adjacencyMatrix)
   inputConfig.adjacencyMatrix = double(~eye(agentNumber)); 
end

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentScenario = scenarioBuilder(agentNumber);
[ agentConfig ] = agentScenario.randomSphere('velocities',inputConfig.velocities,...
                                                 'radius',inputConfig.agentOrbit);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
agentIndex = cell(agentNumber,1);
for index = 1:agentNumber
    agentIndex{index} = inputConfig.agents{index};                      % Get the agents from the input structure
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.quaternion = agentConfig.quaternion(:,index);        % Append properties from the sphereical scenario
    % APPEND THE FORMATION CONTROL ADJACENCY MATRIX
    if isprop(inputConfig.agents{index},'adjacencyMatrix')
        agentIndex{index}.adjacencyMatrix = inputConfig.adjacencyMatrix;
    end
end
                                           
%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
% There are four obstacles in a linear configuration in this example. Their
% placement must be allocated here.

obstacleScenario = scenarioBuilder(obstacleNumber);
% INNER OBSTACLES
[ obstacleConfigA ] = obstacleScenario.randomSphere('velocities',0,...
                                                        'radius',inputConfig.obstacleOrbit);
                                              
                                              
                                              
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};
    nameString = sprintf('OB:%s',obstacleIndex{index}.name);
    obstacleIndex{index}.name = nameString;
    obstacleIndex{index}.VIRTUAL.size = inputConfig.obstacleRadius;        % Size the obstacles 
    % ALLOCATE GLOBAL STATES
    obstacleIndex{index}.globalPosition = obstacleConfigA.position(:,index);
    obstacleIndex{index}.globalVelocity = obstacleConfigA.velocity(:,index);
    obstacleIndex{index}.quaternion = obstacleConfigA.quaternion(:,index);    
end                                                                                             
                                           
                                           
%% ASSIGN GLOBAL PARAMETERS TO THE AGENT INDEX
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
% fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
% agentIndex = cell(objectNumber,1);
% for index = 1:objectNumber
%     agentIndex{index} = scenarioConfig.agents{index};
%     % APPLY GLOBAL STATE VARIABLES
%     agentIndex{index}.globalPosition = agentConfig.position(:,index);
%     agentIndex{index}.globalVelocity = agentConfig.velocity(:,index);
%     agentIndex{index}.quaternion = agentConfig.quaternion(:,index);
% end

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
%     waypointIndex{index}.globalPosition = waypointConfig.position(index,:)';
%     waypointIndex{index}.globalVelocity = waypointConfig.velocity(index,:)';
%     waypointIndex{index}.quaternion = waypointConfig.quaternion(index,:)';
%     waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
% end
% 
% % BUILD THE COLLECTIVE OBJECT INDEX
% objectIndex = horzcat(agentIndex,waypointIndex);

% COLLECT ALL INITIALISED OBSTACLES
objectIndex = vertcat(agentIndex,obstacleIndex);

%% //////////////// CLEAN UP //////////////////////////////////////////////
% SAVE THE FILE
save(inputConfig.file,'objectIndex');
% PLOT THE SCENE
if inputConfig.plot
    agentScenario.plotObjectIndex(objectIndex);                            % Plot the object index
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex