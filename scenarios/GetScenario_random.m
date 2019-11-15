
function [ objectIndex ] = GetScenario_random(varargin)

fprintf('[SCENARIO]\tGetting the random distribution of objects scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct(...
    'file','scenario.mat',...
    'objects',[],...
    'velocity',18,...
    'positionGain',10,...
    'velocityGain',1,...
    'poseGain',pi,...
    'plot',false,....
    'noiseFactor',0);

% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);

% Get the random object configuration
objectConfig = SBinstance.random(...
    'objects',numel(inputConfig.objects),...
    'positionGain',inputConfig.positionGain,...
    'velocityGain',inputConfig.velocityGain,...
    'poseGain',inputConfig.poseGain);

for index = 1:numel(inputConfig.objects)
    % Assign object
    objectIndex{index} = inputConfig.objects{index};
    % APPLY GLOBAL STATE VARIABLES
    objectIndex{index}.SetGLOBAL('position',objectConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1));
    objectIndex{index}.SetGLOBAL('velocity',objectConfig.velocities(:,index) + inputConfig.noiseFactor*randn(3,1));
    objectIndex{index}.SetGLOBAL('quaternion',objectConfig.quaternions(:,index));
end
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end