
function [ objectIndex ] = GetScenario_randomUniform(varargin)

fprintf('[SCENARIO]\tGetting a normally distributed object scenario.\n');

% DEFAULT AGENT/OBSTACLE CONFIGURATION IN THIS EXAMPLE
defaultConfig = struct(...
    'file','scenario.mat',...
    'objects',[],...
    'waypointRadius',0.5,...
    'velocity',18,...
    'positionGain',10,...
    'velocityGain',1,...
    'poseGain',pi,...
    'noiseFactor',0,...
    'is3D',true,...
    'plot',false);

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
    if inputConfig.is3D
        objectIndex{index}.SetGLOBAL('position',objectConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1);
        objectIndex{index}.SetGLOBAL('velocity',objectConfig.velocities(:,index) + inputConfig.noiseFactor*randn(3,1);
        objectIndex{index}.SetGLOBAL('quaternion',objectConfig.quaternions(:,index));
    else
        objectIndex{index}.SetGLOBAL('position',[objectConfig.positions(1:2,index);0] + inputConfig.noiseFactor*[randn(2,1);0];
        objectIndex{index}.SetGLOBAL('velocity',[objectConfig.velocities(1:2,index);0] + inputConfig.noiseFactor*[randn(2,1);0];
        eta = OMAS_geometry.quaternionToEulers(objectConfig.quaternions(:,index));
        objectIndex{index}.SetGLOBAL('quaternion',OMAS_geometry.eulersToQuaternion([eta(3);0;0]));
    end
end
% Get the random object configuration
waypointConfig = SBinstance.randomUniform(...
    'objects',numel(inputConfig.objects),...
    'positionGain',inputConfig.positionGain,...
    'velocityGain',inputConfig.velocityGain,...
    'poseGain',inputConfig.poseGain);

fprintf('[SCENARIO]\tAssigning waypoints to agents.\n'); 
waypointSet = cell(size(inputConfig.objects));
for index = 1:numel(inputConfig.objects)
    % WAYPOINT SET A
    waypointSet{index} = waypoint('radius',inputConfig.waypointRadius,'priority',1,'name',sprintf('WP-%s',objectIndex{index}.name));
    % APPLY GLOBAL STATE VARIABLES
    if inputConfig.is3D
        waypointSet{index}.SetGLOBAL('position',waypointConfig.positions(:,index) + inputConfig.noiseFactor*randn(3,1));
        waypointSet{index}.SetGLOBAL('velocity',waypointConfig.velocities(:,index));
        waypointSet{index}.SetGLOBAL('quaternion',waypointConfig.quaternions(:,index));
    else
        waypointSet{index}.SetGLOBAL('position',[waypointConfig.positions(1:2,index);0] + inputConfig.noiseFactor*[randn(2,1);0]);
        waypointSet{index}.SetGLOBAL('velocity',[waypointConfig.velocities(1:2,index);0] + inputConfig.noiseFactor*[randn(2,1);0]);
        eta = OMAS_geometry.quaternionToEulers(waypointConfig.quaternions(:,index));
        waypointSet{index}.SetGLOBAL('quaternion',OMAS_geometry.eulersToQuaternion([eta(3);0;0]));
    end
    waypointSet{index} = waypointSet{index}.CreateAgentAssociation(objectIndex{index});  % Create waypoint with association to agent
end
% Concatinate the objects and way-points    
objectIndex = [objectIndex;waypointSet];

% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end