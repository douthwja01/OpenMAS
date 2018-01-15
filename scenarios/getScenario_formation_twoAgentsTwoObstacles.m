function [ objectIndex ] = getScenario_formation_twoAgentsTwoObstacles(varargin)
% This function generates the two agent, two obstacle example, termed
% scenario A in the formation control/ collision avoidance study.

fprintf('[SCENARIO]\tGetting the two agent, two obstacle formation control example.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',10,...
                       'obstacles',2,...
                       'obstacleRadius',5,...
                       'obstacleOrbit',10,...
                       'velocities',0,...
                       'adjacencyMatrix',[],...                            % The globally specified adjacency matrix
                       'plot',0);                     
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);

% THIS AGENT REQUIRES ONLY TWO AGENTS, CONFIRM
agentNumber = numel(inputConfig.agents);
assert(agentNumber == 2,'This scenario (scenario A) requires two input agents, specified by the "agent" attribute.');

if isnumeric(inputConfig.obstacles)
    obstacleSet = cell(inputConfig.obstacles,1);
    for index = 1:inputConfig.obstacles
       obstacleSet{index} = obstacle();
    end
    inputConfig.obstacles = obstacleSet;
end

% DECLARE THE NUMBER OF OBSTACLES
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
    agentIndex{index}.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.quaternion = agentConfig.quaternion(:,index);        % Append properties from the sphereical scenario
    % APPEND THE FORMATION CONTROL ADJACENCY MATRIX
    if isprop(inputConfig.agents{index},'adjacencyMatrix')
        agentIndex{index}.adjacencyMatrix = inputConfig.adjacencyMatrix;
    end
end

%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
obstacleScenario = scenarioBuilder(obstacleNumber);
[ obstacleConfig ] = obstacleScenario.planarRing('radius',inputConfig.obstacleOrbit,...
                                               'velocity',0,...
                                            'offsetAngle',pi);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};              % Get the agents from the input structure
    obstacleIndex{index}.name = sprintf('OB:%s',inputConfig.obstacles{index}.name);
    obstacleIndex{index}.VIRTUAL.size = inputConfig.obstacleRadius;
    % APPLY GLOBAL STATE VARIABLES
    obstacleIndex{index}.globalPosition = obstacleConfig.position(:,index);
    obstacleIndex{index}.globalVelocity = obstacleConfig.velocity(:,index);
    obstacleIndex{index}.quaternion = obstacleConfig.quaternion(:,index);  % Append properties from the sphereical scenario
end

%% /////////////// CLEAN UP ///////////////////////////////////////////////

% BUILD THE COMPLETE OBJECT SET
objectIndex = vertcat(agentIndex,obstacleIndex); 
% SAVE THE FILE
save(inputConfig.file,'objectIndex');

% PLOT THE SCENE
if inputConfig.plot
    agentScenario.plotObjectIndex(objectIndex);                            % Plot the object index
end

% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end