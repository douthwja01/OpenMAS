function [ objectIndex ] = GetScenario_formation_fourObstacles(varargin)
% This function generates the four agent, four obstacle example, termed
% scenario D in the formation control/ collision avoidance study.

fprintf('[SCENARIO]\tGetting the four agent, four obstacle formation control example.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentOrbit',10,...
    'agentVelocity',0,...
    'obstacles',4,...
    'obstacleRadius',1,...
    'obstacleOrbit',5,...
    'adjacencyMatrix',[],...                            % The globally specified adjacency matrix
    'plot',false);                     
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);

% AGENT CONDITIONING
agentNumber = numel(inputConfig.agents);
assert(agentNumber == 4,'This scenario requires four input agents, specified by the "agent" attribute.');

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
% INNER AGENTS
[ agentConfigA ] = SBinstance.planarRing(...
    'objects',agentNumber/2,...
    'velocities',inputConfig.agentVelocity,...
    'radius',inputConfig.agentOrbit);
% OUTER AGENTS                                           
[ agentConfigB ] = SBinstance.planarRing(...
    'objects',agentNumber/2,...
    'velocities',inputConfig.agentVelocity,...
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
agentIndex{1}.SetGLOBAL('position',agentConfigB.positions(:,1));
agentIndex{1}.SetGLOBAL('velocity',agentConfigB.velocities(:,1));
agentIndex{1}.SetGLOBAL('quaternion',agentConfigB.quaternions(:,1));          % Append properties from the sphereical scenario
agentIndex{2}.SetGLOBAL('position',agentConfigA.positions(:,1));
agentIndex{2}.SetGLOBAL('velocity',agentConfigA.velocities(:,1));
agentIndex{2}.SetGLOBAL('quaternion',agentConfigA.quaternions(:,1));          % Append properties from the sphereical scenario
agentIndex{3}.SetGLOBAL('position',agentConfigA.positions(:,2));
agentIndex{3}.SetGLOBAL('velocity',agentConfigA.velocities(:,2));
agentIndex{3}.SetGLOBAL('quaternion',agentConfigA.quaternions(:,2)); 
agentIndex{4}.SetGLOBAL('position',agentConfigB.positions(:,2));
agentIndex{4}.SetGLOBAL('velocity',agentConfigB.velocities(:,2));
agentIndex{4}.SetGLOBAL('quaternion',agentConfigB.quaternions(:,2));                                                 

%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
% The four obstacles are positioned in a ring around the center
[ obstacleConfig ] = SBinstance.planarRing(...
    'objects',obstacleNumber,...
    'radius',inputConfig.obstacleOrbit,...
    'offsetAngle',pi,...
    'velocity',0);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(obstacleNumber,1);
for index = 1:obstacleNumber
    obstacleIndex{index} = inputConfig.obstacles{index};                                    % Get the agents from the input structure
    obstacleIndex{index}.name = sprintf('OB-%s',inputConfig.obstacles{index}.name);
    obstacleIndex{index}.radius = inputConfig.obstacleRadius;
    % APPLY GLOBAL STATE VARIABLES
    obstacleIndex{index}.SetGLOBAL('position',obstacleConfig.positions(:,index));
    obstacleIndex{index}.SetGLOBAL('velocity',obstacleConfig.velocities(:,index));
    obstacleIndex{index}.SetGLOBAL('quaternion',obstacleConfig.quaternions(:,index));  % Append properties from the sphereical scenario
end

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = vertcat(agentIndex,obstacleIndex); 
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);                            % Plot the object index
end
clearvars -except objectIndex % Clean-up
end