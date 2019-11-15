function [ objectIndex ] = GetScenario_formation_split(varargin)
% This function generates the four agent, four obstacle example, termed
% scenario D in the formation control/ collision avoidance study.

fprintf('[SCENARIO]\tGetting the four agent, four obstacle formation control example.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct(...
    'file','scenario.mat',...
    'origin',zeros(3,1),...
    'agents',[],...
    'agentVelocity',0,...
    'agentSpacing',5,...
    'waypointRadius',1,...
    'plot',false,...
    'noiseFactor',0,...
    'adjacencyMatrix',[]); 

% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);                       
% Input sanity check
inputConfig.agents = reshape(inputConfig.agents,[numel(inputConfig.agents),1]); % Ensure formatting
assert(~isempty(inputConfig.agents),'Please provide an agent vector.');

% DESIGN THE ADJACENCY MATRIX
if isempty(inputConfig.adjacencyMatrix)
    fprintf('[SCENARIO]\tUsing the default adjacency matrix.');
    % The adjacency matrix describes the desired separation between each
    % agents. 
    inputConfig.adjacencyMatrix = double(~eye(numel(inputConfig.agents)));  % DESIGN AN ADJACENCY MATRIX
end

% /////////////////// DESIGN THE OBSTACLE CONFIGURATION ///////////////////
% This scenario has a constant obstacle configuration irrespect of the
% objects in the field.
obstacleIndex = cell(1);
obstacleIndex{1} = obstacle_cuboid('name','Divider','Xscale',10,'Zscale',10);
obstacleIndex{1}.SetGLOBAL('position',[20;0;0]);

% ///////////////////// DESIGN THE AGENT CONFIGURATION ////////////////////
agentNumber = numel(inputConfig.agents);
agentIndex  = inputConfig.agents;
% Define the agent scenario
diskConfig = SBinstance.planarDisk(...
    'objects',agentNumber,...
    'pointA',inputConfig.origin - [0;0;1],...
    'pointB',inputConfig.origin,...              % Center of the disk
    'scale',inputConfig.agentSpacing);

for index = 1:agentNumber 
    % ASSIGN THE GLOBAL ADJACENCY MATRIX 
    if isprop(agentIndex{index},'adjacencyMatrix')
        agentIndex{index}.adjacencyMatrix = inputConfig.adjacencyMatrix;
    end
    % ASSIGN GLOBAL PROPERTIES
    agentIndex{index}.SetGLOBAL('position',diskConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1));
    agentIndex{index}.SetGLOBAL('velocity',diskConfig.velocities(:,index));
    agentIndex{index}.SetGLOBAL('quaternion',diskConfig.quaternions(:,index));
end

% ////////// DESIGN THE REPRESENTATIVE WAYPOINT CONFIGURATION /////////////
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
waypointIndex = cell(agentNumber,1);
for index = 1:agentNumber
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index,1} = waypoint('radius',inputConfig.waypointRadius,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index,1}.SetGLOBAL('position',[40;0;0]);
    waypointIndex{index,1}.SetGLOBAL('velocity',zeros(3,1));
    waypointIndex{index,1} = waypointIndex{index}.CreateAgentAssociation(agentIndex{index},5);  % Create waypoint with association to agent
end

% /////////////////////////////// CLEAN UP ////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex;obstacleIndex;waypointIndex]; 
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);                          % Plot the object index
end
clearvars -except objectIndex % Clean-ups
end