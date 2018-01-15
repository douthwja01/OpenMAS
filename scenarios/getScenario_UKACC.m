function [ objectIndex ] = getScenario_UKACC(varargin)
% This function generates the evalulation scenario intended for UKACC 2018.
% The scenario defines a generic concentric scenario for a given agent set,
% no obstacles are used as merely inter-agent avoidance. 

% James A. Douthwaite 04/01/2018

fprintf('[SCENARIO]\tGetting UKACC cocentric evaluation scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',15,...
                       'agentRadius',1,...
                       'waypoints',[],...
                       'waypointOrbit',10,...
                       'waypointRadius',0.5,...
                       'velocities',0,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex = inputConfig.agents;                                           % Declare the agent set
agentNumber = numel(agentIndex);                                           % Declare the number of agents

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentScenario = scenarioBuilder(agentNumber);
[agentConfig] = agentScenario.planarRing('velocities',inputConfig.velocities,...
                                             'radius',inputConfig.agentOrbit);   

%% REBUILD THE AGENT INDEX UNDER THIS CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent definitions...\n'); 
for index = 1:agentNumber
    % MODIFY AGENT SIZE
    agentIndex{index}.VIRTUAL.size = inputConfig.agentRadius;
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.quaternion = agentConfig.quaternion(:,index);
end

%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% DEFINE THE WAYPOINT CONFIGURATIONS
waypointPlanarRotation = pi;
[ waypointConfig] = agentScenario.planarRing('radius',inputConfig.waypointOrbit,...
                                         'velocities',0,...
                                          'zeroAngle',waypointPlanarRotation);
                                        
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
waypointIndex = cell(size(agentIndex));
for index = 1:agentNumber
    nameString = sprintf('WP:%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.globalPosition = waypointConfig.position(:,index);
    waypointIndex{index}.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.quaternion = waypointConfig.quaternion(:,index);
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{index});  % Create waypoint with association to agent
end

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentIndex,waypointIndex]; 
% SAVE THE FILE
save(inputConfig.file,'objectIndex');
% PLOT THE SCENE
if inputConfig.plot
    agentScenario.plotObjectIndex(objectIndex);                            % Plot the object index
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end