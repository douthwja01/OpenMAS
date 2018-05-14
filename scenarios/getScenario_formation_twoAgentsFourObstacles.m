function [ objectIndex ] = getScenario_formation_twoAgentsFourObstacles(varargin)
% This function generates the two agent, four obstacle example, termed
% scenario B in the formation control/ collision avoidance study.

fprintf('[SCENARIO]\tGetting the two agent, four obstacle formation control example.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',20,...
                       'obstacles',4,...
                       'obstacleRadius',2,...
                       'obstacleOrbit',10,...
                       'velocities',0,...
                       'adjacencyMatrix',[],...                            % The globally specified adjacency matrix
                       'plot',0);                     
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);

% THIS AGENT REQUIRES ONLY TWO AGENTS, CONFIRM
agentNumber = numel(inputConfig.agents);
assert(agentNumber == 2,'This scenario (scenario B) requires two input agents, specified by the "agent" attribute.');
if isnumeric(inputConfig.obstacles)
    obstacleSet = cell(inputConfig.obstacles,1);
    for index = 1:inputConfig.obstacles
       obstacleSet{index} = obstacle();
    end
    inputConfig.obstacles = obstacleSet;
end
obstacleNumber = numel(inputConfig.obstacles);  

% DESIGN THE DESIRED SEPERATION MATRIX (ADJACENCY MATRIX)
% The adjacency matrix is indexed by objectID in scenarioConfig.adjacencyMatrix
% OTHERWISE ASSIGN DEFAULT ADJACENCY MATRIX
if isempty(inputConfig.adjacencyMatrix)
   inputConfig.adjacencyMatrix = double(~eye(agentNumber)); 
end

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentScenario = scenarioBuilder(agentNumber);
[ agentConfig ] = agentScenario.planarRing('velocities',inputConfig.velocities,...
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
spacingCoefficient = 2;
obstacleScenario = scenarioBuilder(obstacleNumber/2);
% INNER OBSTACLES
[ obstacleConfigA ] = obstacleScenario.planarRing('velocities',0,...
                                                   'radius',inputConfig.obstacleOrbit);
% OUTER AGENTS                                           
[ obstacleConfigB ] = obstacleScenario.planarRing('velocities',0,...
                                                   'radius',inputConfig.obstacleOrbit*spacingCoefficient);
                                               
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};
    nameString = sprintf('OB:%s',obstacleIndex{index}.name);    
    obstacleIndex{index}.name = nameString;
    obstacleIndex{index}.VIRTUAL.size = inputConfig.obstacleRadius;        % Size the obstacles 
end                                                   
% ASSIGN AGENT GLOBAL PROPERTIES, ONE SIDE OF THE RINGS TO THE OTHER
obstacleIndex{1}.VIRTUAL.globalPosition = obstacleConfigB.position(:,1);
obstacleIndex{1}.VIRTUAL.globalVelocity = obstacleConfigB.velocity(:,1);
obstacleIndex{1}.VIRTUAL.quaternion = obstacleConfigB.quaternion(:,1);             % Append properties from the sphereical scenario
obstacleIndex{2}.VIRTUAL.globalPosition = obstacleConfigA.position(:,1);
obstacleIndex{2}.VIRTUAL.globalVelocity = obstacleConfigA.velocity(:,1);
obstacleIndex{2}.VIRTUAL.quaternion = obstacleConfigA.quaternion(:,1);             % Append properties from the sphereical scenario
obstacleIndex{3}.VIRTUAL.globalPosition = obstacleConfigA.position(:,2);
obstacleIndex{3}.VIRTUAL.globalVelocity = obstacleConfigA.velocity(:,2);
obstacleIndex{3}.VIRTUAL.quaternion = obstacleConfigA.quaternion(:,2); 
obstacleIndex{4}.VIRTUAL.globalPosition = obstacleConfigB.position(:,2);
obstacleIndex{4}.VIRTUAL.globalVelocity = obstacleConfigB.velocity(:,2);
obstacleIndex{4}.VIRTUAL.quaternion = obstacleConfigB.quaternion(:,2);                                                  

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
end