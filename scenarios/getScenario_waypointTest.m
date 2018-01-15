function [ objectIndex ] = getScenario_waypointTest(varargin)
% This scenario is designed to present a waypoint following exercise.

fprintf('[SCENARIO]\tGetting the waypoint following exercise.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario_waypointTest.mat',...
                       'agents',[],...
                       'noiseFactor',0,...
                       'plot',0);  
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);
agentIndex = scenarioConfig.agents;
agentNumber = numel(agentIndex);                    % Declare the number of agents

%% DEFINE THE AGENT CONFIGURATION
% MOVE THROUGH THE AGENTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning agent definition...\n'); 
for i=1:agentNumber
    agentIndex{1}.globalPosition = [0;0;0] + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{1}.globalVelocity = [10;0;0] + scenarioConfig.noiseFactor*randn(3,1);
    agentIndex{1}.quaternion = [1;0;0;0];
end

%% DEFINE THE WAYPOINT CONFIGURATION
fprintf('[SCENARIO]\tBuilding the new scenario...\n');
waypointNumber = 3;
testScenario = scenarioBuilder(waypointNumber);
angleOffset = -(10/9)*pi;
[ waypointConfig] = testScenario.planarAngle('radius',200,...
                                             'pointA',[200;-200;-1],...
                                             'pointB',[200;-200;0],...
                                             'velocities',0,...
                                             'zeroAngle',angleOffset);

% MOVE THROUGH THE WAYPOINTS AND INITIALISE WITH GLOBAL PROPERTIES
fprintf('[SCENARIO]\tAssigning waypoint definitions:\n'); 
for index = 1:waypointNumber
    waypointIndex{index} = waypoint('radius',2);
    % APPLY GLOBAL STATE VARIABLES
    waypointIndex{index}.globalPosition = waypointConfig.position(index,:)';
    waypointIndex{index}.globalVelocity = waypointConfig.velocity(index,:)';
    waypointIndex{index}.quaternion = waypointConfig.quaternion(index,:)';
    waypointPriority = 1/index;
    waypointIndex{index} = waypointIndex{index}.createAgentAssociation(agentIndex{1},waypointPriority);  % Create waypoint with association to agent
end

% BUILD THE COLLECTIVE OBJECT INDEX
objectIndex = horzcat(agentIndex,waypointIndex);
% PLOT THE SCENE
if scenarioConfig.plot
    scenarioBuilder.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex','agentIndex','waypointIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex

end

% PARSE CONFIGURATION INPUTS
function [agentIndex,noiseFactor,configFile] = parseInputs(inputSet)

% DEFAULT FILE NAME
configFile = 'scenario.mat';
noiseFactor = 0;
agentIndex = {};

% VISABILTY REDUCTION RATIO (IN FUTURE PROPORTIONAL TO ISO PROPERTIES)
tmp = strncmpi(inputSet,'filename',4);
if any(tmp)
    configFile = inputSet{find(tmp) + 1};
end

tmp = strncmpi(inputSet,'agents',5)|strncmpi(inputSet,'objects',6); 
if any(tmp)
    agentIndex = inputSet{find(tmp) + 1}; 
    assert(iscell(agentIndex) == 1,'Figure set must be specified as a cell array of string names.');
end

tmp = strncmpi(inputSet,'noisefactor',5); % |strncmpi(inputSet,'objects',6); 
if any(tmp)
    noiseFactor = inputSet{find(tmp) + 1}; 
    assert(isnumeric(noiseFactor) == 1,'Noise factor must be a numeric factor.');
end

% CHECK THERE IS SOMETHING TO SIMULATE
if isempty(agentIndex)
    error('Please provide a cell array of agent objects to simulation');
end

end