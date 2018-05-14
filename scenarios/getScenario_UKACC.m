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
                       'agentVelocity',0,...
                       'offsetAngle',0,...
                       'agentRadius',0.5,...    % Width of 1m
                       'waypoints',[],...
                       'waypointOrbit',15,...
                       'waypointOffsetAngle',[],...
                       'waypointRadius',0.5,... % Diameter of 1m
                       'noiseFactor',0,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex = inputConfig.agents;                                           % Declare the agent set
agentNumber = numel(agentIndex);                                           % Declare the number of agents

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
agentScenario = scenarioBuilder(agentNumber);
[agentConfig] = agentScenario.planarRing('velocities',inputConfig.agentVelocity,...
                                             'radius',inputConfig.agentOrbit);   

%% REBUILD THE AGENT INDEX UNDER THIS CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:agentNumber
    % MODIFY AGENT SIZE
    agentIndex{index}.VIRTUAL.radius = inputConfig.agentRadius;
    % APPLY GLOBAL STATE VARIABLES
    agentIndex{index}.VIRTUAL.globalPosition = agentConfig.position(:,index);
    agentIndex{index}.VIRTUAL.globalVelocity = agentConfig.velocity(:,index);
    agentIndex{index}.VIRTUAL.quaternion = agentConfig.quaternion(:,index);
end

%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% DEFINE THE WAYPOINT CONFIGURATIONS
waypointPlanarRotation = pi;
[ waypointConfig] = agentScenario.planarRing('velocities',0,...
                                                 'radius',inputConfig.waypointOrbit,...
                                              'zeroAngle',waypointPlanarRotation);
                                        
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoints to agents.\n'); 
waypointIndex = cell(size(agentIndex));
for index = 1:agentNumber
    nameString = sprintf('WP-%s',agentIndex{index}.name);
    waypointIndex{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
%     waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(:,index) + inputConfig.noiseFactor*randn(3,1);
    waypointIndex{index}.VIRTUAL.globalPosition = waypointConfig.position(:,index) + inputConfig.noiseFactor*randn(3,1);
    waypointIndex{index}.VIRTUAL.globalVelocity = waypointConfig.velocity(:,index);
    waypointIndex{index}.VIRTUAL.quaternion = waypointConfig.quaternion(:,index);
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
% 
fprintf('[SCENARIO]\tDone.\n'); 
end