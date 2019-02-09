function [ objectIndex ] = getScenario_earthOrbit(varargin)
% This function is designed to generate a scenario that represents a
% satellite orbiting earth.

fprintf('[SCENARIO]\tGetting the Earth and satillite scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct('file','scenario.mat',...
                       'agents',[],...
                       'agentOrbit',435E3,...       % The true orbit of the ISS
                       'agentVelocity',0,...
                       'scale',5E3,...
                       'waypointOrbit',[],...
                       'waypointOffsetAngle',[],...
                       'waypointRadius',0.5,...
                       'noiseFactor',0,...
                       'plot',0);  
                   
% PARSE THE USER OVERRIDES USING THE SCENARIO BUILDER
[scenarioConfig] = scenarioBuilder.configurationParser(defaultConfig,varargin);

% GENERATE THE EARTH AS THE FOCAL POINT
objectIndex{1} = earth('name','Earth'); % Will be initiallised at the origin.

% ////////////////// DEFINE THE MOONS'S ORBITAL PROPERTIES ////////////////
objectIndex{2} = moon('name','Moon');
orbitalInclination = deg2rad(objectIndex{2}.inclination);
axisProjection = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);
% GENERATE THE MOON IN ORBIT
ISS_scenario = scenarioBuilder(1);
Moon_config = ISS_scenario.planarRing('objects',1,'radius',objectIndex{2}.orbit,...
                                      'pointA',axisProjection,'pointB',zeros(3,1));
objectIndex{2}.VIRTUAL.globalVelocity = [0;objectIndex{2}.orbitalSpeed;0];                    % Initialise with tangential oribit speed of 7.67km/s
objectIndex{2}.VIRTUAL.globalPosition = Moon_config.position;               % Initialise at the apogee (403k), perigee is 406km
objectIndex{2}.VIRTUAL.quaternion = Moon_config.quaternion;

% ////////////////// DEFINE THE ISS'S ORBITAL PROPERTIES //////////////////
% GENERATE THE ISS IN ORBIT
objectIndex{3} = ISS('name','ISS','scale',scenarioConfig.scale);
orbitalInclination = deg2rad(0);
axisProjection = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);
ISS_scenario = scenarioBuilder(1);
ISS_config = ISS_scenario.planarRing('objects',1,'radius',objectIndex{3}.orbit,...
                                      'pointA',axisProjection,'pointB',zeros(3,1));
objectIndex{3}.VIRTUAL.globalVelocity = [0;objectIndex{3}.orbitalSpeed;0];                      % Initialise with tangential oribit speed of 7.67km/s
objectIndex{3}.VIRTUAL.globalPosition = ISS_config.position;               % Initialise at the apogee (403k), perigee is 406km
objectIndex{3}.VIRTUAL.quaternion = ISS_config.quaternion;

% PLOT THE SCENE
if scenarioConfig.plot
    ISS_scenario.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end