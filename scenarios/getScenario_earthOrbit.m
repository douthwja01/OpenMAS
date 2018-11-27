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

% DEFINE THE ISS'S ORBITAL PROPERTIES
orbitalAltitude = 6.371E+06 + 406E3; 
orbitalInclination = deg2rad(51.64);
axisProjection = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);

ISS_scenario = scenarioBuilder(1);
ISS_config = ISS_scenario.planarRing('objects',1,...
                                      'radius',orbitalAltitude,...
                                      'pointA',axisProjection,...
                                      'pointB',zeros(3,1));

% GENERATE THE ISS IN ORBIT
objectIndex{2} = ISS('name','ISS','scale',scenarioConfig.scale);
objectIndex{2}.VIRTUAL.globalVelocity = [0;7.67E3;0];                      % Initialise with tangential oribit speed of 7.67km/s
objectIndex{2}.VIRTUAL.globalPosition = ISS_config.position;               % Initialise at the apogee (403k), perigee is 406km
objectIndex{2}.VIRTUAL.quaternion = ISS_config.quaternion;



% PLOT THE SCENE
if scenarioConfig.plot
    ISS_scenario.plotObjectIndex(objectIndex);
end
% SAVE THE FILE
save(scenarioConfig.file,'objectIndex');
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end