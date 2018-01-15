function [ objectIndex ] = getScenario_formation_fourAgentsFourObstacles(varargin)
% This function generates the four agent, four obstacle example, termed
% scenario D in the formation control/ collision avoidance study.

fprintf('[SCENARIO]\tGetting the four agent, four obstacle formation control example.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',10,...
                       'obstacles',4,...
                       'obstacleRadius',2,...
                       'obstacleOrbit',10,...
                       'velocities',0,...
                       'adjacencyMatrix',[],...                            % The globally specified adjacency matrix
                       'plot',0);                     
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);

% AGENT CONDITIONING
agentNumber = numel(inputConfig.agents);
assert(agentNumber == 4,'This scenario (scenario D) requires four input agents, specified by the "agent" attribute.');

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
agentScenario = scenarioBuilder(agentNumber/2);
% INNER AGENTS
[ agentConfigA ] = agentScenario.planarRing('velocities',inputConfig.velocities,...
                                                'radius',inputConfig.agentOrbit);
% OUTER AGENTS                                           
[ agentConfigB ] = agentScenario.planarRing('velocities',inputConfig.velocities,...
                                                'radius',inputConfig.agentOrbit*1.5);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
agentIndex = cell(agentNumber,1);
for index = 1:agentNumber
    agentIndex{index} = inputConfig.agents{index};
    % APPEND THE FORMATION CONTROL ADJACENCY MATRIX
    if isprop(inputConfig.agents{index},'adjacencyMatrix')
        agentIndex{index}.adjacencyMatrix = inputConfig.adjacencyMatrix;
    end
end
% ASSIGN AGENT GLOBAL PROPERTIES, ONE SIDE OF THE RINGS TO THE OTHER
agentIndex{1}.globalPosition = agentConfigB.position(:,1);
agentIndex{1}.globalVelocity = agentConfigB.velocity(:,1);
agentIndex{1}.quaternion = agentConfigB.quaternion(:,1);                   % Append properties from the sphereical scenario
agentIndex{2}.globalPosition = agentConfigA.position(:,1);
agentIndex{2}.globalVelocity = agentConfigA.velocity(:,1);
agentIndex{2}.quaternion = agentConfigA.quaternion(:,1);                   % Append properties from the sphereical scenario
agentIndex{3}.globalPosition = agentConfigA.position(:,2);
agentIndex{3}.globalVelocity = agentConfigA.velocity(:,2);
agentIndex{3}.quaternion = agentConfigA.quaternion(:,2); 
agentIndex{4}.globalPosition = agentConfigB.position(:,2);
agentIndex{4}.globalVelocity = agentConfigB.velocity(:,2);
agentIndex{4}.quaternion = agentConfigB.quaternion(:,2);                                                 

%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
% The four obstacles are positioned in a ring around the center
obstacleScenario = scenarioBuilder(obstacleNumber);
[ obstacleConfig ] = obstacleScenario.planarRing('radius',inputConfig.obstacleOrbit,...
                                            'offsetAngle',pi,...
                                               'velocity',0);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};                                    % Get the agents from the input structure
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