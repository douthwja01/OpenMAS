function [ objectIndex ] = getScenario_fourCuboidObstacles(varargin)
% This function is designed to generate a simple cube of obstacles such
% that the agents must navigate around them simultaneously in order to
% reach their desired goals.

fprintf('[SCENARIO]\tGetting Simple thre obstacle, one waypoint scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...   
                       'agentVelocity',0,...
                       'agentRadius',0.5,...    % Width of 1m
                       'agentOrbit',10,...
                       'obstacleRadius',1,...
                       'obstacleOrbit',5,...
                       'waypoints',[],...
                       'waypointRadius',0.5,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
% HANDLE OBSTACLE GENERATION/OBSTACLES PROVIDED

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentScenario = scenarioBuilder(numel(inputConfig.agents));
% INNER AGENTS
[ agentConfig ] = agentScenario.planarRing('velocities',inputConfig.agentVelocity,...
                                                'radius',inputConfig.agentOrbit);

% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
agentIndex = cell(1,numel(inputConfig.agents));
for index = 1:numel(inputConfig.agents)
    agentIndex{index} = inputConfig.agents{index};                                    % Get the agents from the input structure
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);  % Append properties from the sphereical scenario
end  
 
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint global parameters...\n'); 
waypointIndex = cell(1,numel(inputConfig.agents));
for index = 1:numel(inputConfig.agents)
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.VIRTUAL.globalPosition =  1.2*agentConfig.position(:,index);
    waypointIndex{index}.VIRTUAL.globalVelocity = -agentConfig.velocity(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);  % Append properties from the sphereical scenario
    % GENERATE WAY-POINT ASSOCIATION
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index});  % Create waypoint with association to agent
end 

%% //////////////// BUILD THE OBSTACLES GLOBAL STATES /////////////////////
% GET THE SCENARIO BUILDING TOOLS
obstacleScenario = scenarioBuilder(4);
% DEFINE THE AGENT CONFIGURATIONS
[ obstacleConfig ] = obstacleScenario.planarRing('radius',inputConfig.obstacleOrbit,...
                                                 'zeroAngle',0);
                                         
% MOVE THROUGH THE OBSTACLES AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning obstacle global parameters...\n'); 
obstacleIndex = cell(1,obstacleConfig.objects);
for index = 1:obstacleConfig.objects
    label = sprintf('OB-%.0f',index);
    obstacleIndex{index} = obstacle_cuboid('name',label);                                    % Get the agents from the input structure

    % APPLY GLOBAL STATE VARIABLES
    obstacleIndex{index}.VIRTUAL.globalPosition = obstacleConfig.position(:,index);
    obstacleIndex{index}.VIRTUAL.globalVelocity = obstacleConfig.velocity(:,index);
    obstacleIndex{index}.VIRTUAL.quaternion = obstacleConfig.quaternion(:,index);  % Append properties from the sphereical scenario
end                                       


%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex,obstacleIndex,waypointIndex]; 
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


function [obstacleVerticies] = generateObstacleCuboid(cubeCentroid,cubeScale,obstacleRadii)

cubeMapping = 0.5*[-1 -1 -1;
                    1 -1 -1;
                    1  1 -1;
                   -1  1 -1;
                    1  1  1;
                   -1  1  1;
                   -1 -1  1;
                    1 -1  1]'*cubeScale + cubeCentroid; % Relative square positions
                
% GENERATE OBSTACLE CUBOID
obstacleVerticies = cell(1,size(cubeMapping,1));
for vert = 1:size(cubeMapping,2)
    point = cubeMapping(:,vert);
    % GENERATE OBSTACLE
    obstacleVerticies{vert} = obstacle();
    % APPLY GLOBAL STATE VARIABLES
    obstacleVerticies{vert}.VIRTUAL.globalPosition = point';
    obstacleVerticies{vert}.VIRTUAL.globalVelocity = [0;0;0];
    obstacleVerticies{vert}.VIRTUAL.quaternion = [1;0;0;0];  % Append properties from the sphereical scenario
    % APPLY GENERIC GLOBAL STATE VARIABLES
%     obstacleIndex{vert}.name = sprintf('OB-%s',inputConfig.obstacles{index}.name);
    obstacleVerticies{vert}.VIRTUAL.radius = obstacleRadii;
end    

end