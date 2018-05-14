function [ objectIndex ] = getScenario_formation_randomAgentsRandomObstacles(varargin)
% This function designs a scenario where the agents are randomly
% posistioned around a outer sphere, where the obstacles are also
% positioned randomly around a smaller inner sphere.

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
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);        % Append properties from the sphereical scenario
    % APPEND THE FORMATION CONTROL ADJACENCY MATRIX
    if isprop(inputConfig.agents{index},'adjacencyMatrix')
        agentIndex{index}.adjacencyMatrix = inputConfig.adjacencyMatrix;
    end
end
                                           
%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
% There are four obstacles in a linear configuration in this example. Their
% placement must be allocated here.
obstacleScenario = scenarioBuilder(obstacleNumber);
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
    obstacleIndex{index}.VIRTUAL.globalPosition = obstacleConfigA.position(:,index);
    obstacleIndex{index}.VIRTUAL.globalVelocity = obstacleConfigA.velocity(:,index);
    obstacleIndex{index}.VIRTUAL.quaternion = obstacleConfigA.quaternion(:,index);    
end                                                                                             
                                           
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