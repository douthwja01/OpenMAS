function [ objectIndex ] = getScenario_US_twoLines(varargin)
% This scenario constructs the second evaluation scenario 

fprintf('[SCENARIO]\tGetting US two-lines scenario.\n');

%% SCENARIO INPUT HANDLING ////////////////////////////////////////////////
% DEFAULT INPUT CONDITIONS
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'separation',10,...   
                       'agentVelocity',0,...
                       'agentRadius',0.5,...    % Diameter of 1m
                       'waypoints',[],...
                       'waypointRadius',0.5,... % Diameter of 1m
                       'noiseFactor',0,...
                       'plot',0);
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[inputConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
% AGENT CONDITIONS
agentIndex = inputConfig.agents;                                           % Declare the agent set
agentNumber = numel(agentIndex);                                           % Declare the number of agents

% WE DIVID THE NUMBER OF AGENTS INTO TWO GROUPS, ONE FOR EACH LINE, AND
% CREATE ASSOCIATED PAIRS OF WAYPOINTS
assert(mod(agentNumber,2) == 0,'There must be an equal number of agents.');

% SPLIT THE AGENTS INTO GROUPS
agentSetA = agentIndex(1:agentNumber/2);
agentSetB = agentIndex((agentNumber/2)+1:end);

% GLOBAL CONFIGURATION PARAMETERS
paddingMultiplier = 3;
xSpacing = paddingMultiplier*inputConfig.agentRadius;
xLimit = (agentNumber/2)*xSpacing;
agentLineSeparation = 10;
waypointLintSeparation = 12;
setAHeadings =  [0;1;0];
setBHeadings = -[0;1;0];

%% /////////////////// BUILD THE AGENTS GLOBAL STATES /////////////////////
% SCENARIO BUILDER OBJECT 
agentScenario = scenarioBuilder((agentNumber/2));
% CONFIGURATION FOR THE FIRST SET OF AGENTS
[agentConfigA] = agentScenario.line('objects',agentSetA,...
                                    'pointA',[-xLimit;-agentLineSeparation;0],...
                                    'pointB',[ xLimit;-agentLineSeparation;0],...
                                   'heading',setAHeadings,...
                                'velocities',inputConfig.agentVelocity,...
                                    'radius',inputConfig.agentRadius); 
[agentConfigB] = agentScenario.line('objects',agentSetB,...
                                    'pointA',[-xLimit;agentLineSeparation;0],...
                                    'pointB',[ xLimit;agentLineSeparation;0],...
                                   'heading',setBHeadings,...
                                'velocities',inputConfig.agentVelocity,...
                                    'radius',inputConfig.agentRadius); 

%% REBUILD THE AGENT INDEX UNDER THIS CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent global parameters...\n'); 
for index = 1:(agentNumber/2)
    % AGENT SET A
    agentSetA{index}.VIRTUAL.radius = inputConfig.agentRadius;                % Regulate agent radius
    agentSetA{index}.VIRTUAL.globalPosition = agentConfigA.position(:,index); % APPLY GLOBAL STATE VARIABLES
    agentSetA{index}.VIRTUAL.globalVelocity = agentConfigA.velocity(:,index);
    agentSetA{index}.VIRTUAL.quaternion = agentConfigA.quaternion(:,index);
    % AGENT SET B
    agentSetB{index}.VIRTUAL.radius = inputConfig.agentRadius;                % Regulate agent radius
    agentSetB{index}.VIRTUAL.globalPosition = agentConfigB.position(:,index); % APPLY GLOBAL STATE VARIABLES
    agentSetB{index}.VIRTUAL.globalVelocity = agentConfigB.velocity(:,index);
    agentSetB{index}.VIRTUAL.quaternion = agentConfigB.quaternion(:,index);
end

%% //////////////// BUILD THE WAYPOINT GLOBAL STATES //////////////////////
% DEFINE THE WAYPOINT CONFIGURATIONS
waypointScenario = scenarioBuilder((agentNumber/2));
% CONFIGURATION FOR THE FIRST SET OF OBJECTS
[waypointConfigA] = waypointScenario.line('pointA',[-xLimit;waypointLintSeparation;0],...
                                          'pointB',[xLimit;waypointLintSeparation;0],...
                                          'heading',-setAHeadings,...
                                      'velocities',0,...
                                          'radius',inputConfig.waypointRadius);                                      
[waypointConfigB] = waypointScenario.line('pointA',[-xLimit;-waypointLintSeparation;0],...
                                          'pointB',[xLimit;-waypointLintSeparation;0],...
                                         'heading',-setBHeadings,...
                                      'velocities',0,...
                                          'radius',inputConfig.waypointRadius);    
                                      
                                      
% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoints to agents.\n'); 
waypointSetA = cell(size(agentSetA));
waypointSetB = cell(size(agentSetB));
for index = 1:(agentNumber/2)
    % WAYPOINT SET A
    nameString = sprintf('WP-%s',agentSetA{index}.name);
    waypointSetA{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointSetA{index}.VIRTUAL.globalPosition = waypointConfigA.position(:,index) + inputConfig.noiseFactor*[randn(2,1);0];
    waypointSetA{index}.VIRTUAL.globalVelocity = waypointConfigA.velocity(:,index);
    waypointSetA{index}.VIRTUAL.quaternion = waypointConfigA.quaternion(:,index);
    waypointSetA{index} = waypointSetA{index}.createAgentAssociation(agentSetA{index});  % Create waypoint with association to agent
    % WAYPOINT SET B
    nameString = sprintf('WP-%s',agentSetB{index}.name);
    waypointSetB{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',nameString);
    % APPLY GLOBAL STATE VARIABLES
    waypointSetB{index}.VIRTUAL.globalPosition = waypointConfigB.position(:,index) + inputConfig.noiseFactor*[randn(2,1);0];
    waypointSetB{index}.VIRTUAL.globalVelocity = waypointConfigB.velocity(:,index);
    waypointSetB{index}.VIRTUAL.quaternion = waypointConfigB.quaternion(:,index);
    waypointSetB{index} = waypointSetB{index}.createAgentAssociation(agentSetB{index});  % Create waypoint with association to agent
end

%% /////////////// CLEAN UP ///////////////////////////////////////////////
% BUILD THE COMPLETE OBJECT SET
objectIndex = [agentSetA,agentSetB,waypointSetA,waypointSetB]; 
% SAVE THE FILE
save(inputConfig.file,'objectIndex');
% PLOT THE SCENE
if inputConfig.plot
    waypointScenario.plotObjectIndex(objectIndex);                            % Plot the object index
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
% 
fprintf('[SCENARIO]\tDone.\n'); 
end