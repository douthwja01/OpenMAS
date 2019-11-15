function [ objectIndex ] = GetScenario_twoLines(varargin)
% This scenario constructs the second evaluation scenario 

fprintf('[SCENARIO]\tGetting a scenario defined by two opposing lines of agents.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentSeparation',10,...   
    'agentVelocity',0,...
    'agentRadius',0.5,...    % Diameter of 1m
    'waypoints',[],...
    'waypointSeparation',12,...
    'waypointRadius',0.5,... % Diameter of 1m
    'padding',2,...
    'noiseFactor',0,...
    'plot',false);
                   
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex  = inputConfig.agents;                                           % Declare the agent set
agentNumber = numel(agentIndex);                                           % Declare the number of agents

% WE DIVID THE NUMBER OF AGENTS INTO TWO GROUPS, ONE FOR EACH LINE, AND
% CREATE ASSOCIATED PAIRS OF WAYPOINTS
assert(mod(agentNumber,2) == 0,'There must be an equal number of agents in this scenario.');

% SPLIT THE AGENTS INTO GROUPS
agentSetA = agentIndex(1:agentNumber/2);
agentSetB = agentIndex((agentNumber/2)+1:end);

% GLOBAL CONFIGURATION PARAMETERS
xSpacing = inputConfig.padding; %*inputConfig.agentRadius;
xCoords  = linspace(-(agentNumber/4),(agentNumber/4),(agentNumber/2))*xSpacing;

setAHeadings =  [0;1;0];
setBHeadings = -[0;1;0];

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
% CONFIGURATION FOR THE FIRST SET OF AGENTS
[agentConfigA] = SBinstance.line(...
    'objects',agentNumber/2,...
    'pointA',[xCoords(1);  -inputConfig.agentSeparation;0],...
    'pointB',[xCoords(end);-inputConfig.agentSeparation;0],...
    'heading',setAHeadings,...
    'velocities',inputConfig.agentVelocity,...
    'radius',inputConfig.agentRadius);
[agentConfigB] = SBinstance.line(...
    'objects',agentNumber/2,...
    'pointA',[xCoords(end);inputConfig.agentSeparation;0],...
    'pointB',[xCoords(1);inputConfig.agentSeparation;0],...
    'heading',setBHeadings,...
    'velocities',inputConfig.agentVelocity,...
    'radius',inputConfig.agentRadius);

%% REBUILD THE AGENT INDEX UNDER THIS CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:(agentNumber/2)
    % AGENT SET A
    agentSetA{index}.radius = inputConfig.agentRadius;             % Regulate agent radius
    agentSetA{index}.SetGLOBAL('position',agentConfigA.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]); 
    agentSetA{index}.SetGLOBAL('velocity',agentConfigA.velocities(:,index));
    agentSetA{index}.SetGLOBAL('quaternion',agentConfigA.quaternions(:,index));
    % AGENT SET B
    agentSetB{index}.radius = inputConfig.agentRadius;             % Regulate agent radius
    agentSetB{index}.SetGLOBAL('position',agentConfigB.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]); 
    agentSetB{index}.SetGLOBAL('velocity',agentConfigB.velocities(:,index));
    agentSetB{index}.SetGLOBAL('quaternion',agentConfigB.quaternions(:,index));
end

%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% CONFIGURATION FOR THE FIRST SET OF OBJECTS
[waypointConfigA] = SBinstance.line(...
    'objects',agentNumber/2,...
    'pointA',[xCoords(end);  inputConfig.waypointSeparation;0],...
    'pointB',[xCoords(1);inputConfig.waypointSeparation;0],...
    'heading',-setAHeadings,...
    'velocities',0,...
    'radius',inputConfig.waypointRadius);
[waypointConfigB] = SBinstance.line(...
    'objects',agentNumber/2,...
    'pointA',[xCoords(1);-inputConfig.waypointSeparation;0],...
    'pointB',[xCoords(end);  -inputConfig.waypointSeparation;0],...
    'heading',-setBHeadings,...
    'velocities',0,...
    'radius',inputConfig.waypointRadius);
                                      
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoints to agents.\n'); 
waypointSetA = cell(size(agentSetA)); waypointSetB = cell(size(agentSetB));
for index = 1:(agentNumber/2)
    % WAYPOINT SET A
    waypointSetA{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',sprintf('WP-%s',agentSetA{index}.name));
    % APPLY GLOBAL STATE VARIABLES
    waypointSetA{index}.SetGLOBAL('position',waypointConfigA.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]);
    waypointSetA{index}.SetGLOBAL('velocity',waypointConfigA.velocities(:,index));
    waypointSetA{index}.SetGLOBAL('quaternion',waypointConfigA.quaternions(:,index));
    waypointSetA{index} = waypointSetA{index}.CreateAgentAssociation(agentSetA{index});  % Create waypoint with association to agent
    % WAYPOINT SET B
    waypointSetB{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',sprintf('WP-%s',agentSetB{index}.name));
    % APPLY GLOBAL STATE VARIABLES
    waypointSetB{index}.SetGLOBAL('position',waypointConfigB.positions(:,index) + inputConfig.noiseFactor*[randn(2,1);0]);
    waypointSetB{index}.SetGLOBAL('velocity',waypointConfigB.velocities(:,index));
    waypointSetB{index}.SetGLOBAL('quaternion',waypointConfigB.quaternions(:,index));
    waypointSetB{index} = waypointSetB{index}.CreateAgentAssociation(agentSetB{index});  % Create waypoint with association to agent
end

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentSetA,agentSetB,waypointSetA,waypointSetB]; 
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);                            % Plot the object index
end
clearvars -except objectIndex   % Clean-up
fprintf('[SCENARIO]\tDone.\n'); 
end