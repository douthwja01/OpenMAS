
function [ objectIndex ] = GetScenario_randomUniform(varargin)

fprintf('[SCENARIO]\tGetting a normally distributed object scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct('file','scenario.mat',...
                       'objects',[],...
                       'velocity',18,...
                       'plot',0,...
                       'positionGain',10,...
                       'velocityGain',1,...
                       'poseGain',pi,...
                       'noiseFactor',0);
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);

% Get the random object configuration
objectConfig = SBinstance.randomUniform(...
    'objects',numel(inputConfig.objects),...
    'positionGain',inputConfig.positionGain,...
    'velocityGain',inputConfig.velocityGain,...
    'poseGain',inputConfig.poseGain);

for index = 1:numel(inputConfig.objects)
    % Assign object
    objectIndex{index} = inputConfig.objects{index};
    % APPLY GLOBAL STATE VARIABLES
    objectIndex{index}.VIRTUAL.globalPosition = objectConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1);
    objectIndex{index}.VIRTUAL.globalVelocity = objectConfig.velocities(:,index) + inputConfig.noiseFactor*randn(3,1);
    objectIndex{index}.VIRTUAL.quaternion = objectConfig.quaternions(:,index);
end
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end