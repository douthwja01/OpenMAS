function [ objectIndex ] = GetScenario_earthOrbit(varargin)
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
                   
% Instanciate the scenario builder
SBinstance = scenarioBuilder();
% Parse user inputs 
[inputConfig] = SBinstance.configurationParser(defaultConfig,varargin);  

% GENERATE THE EARTH AS THE FOCAL POINT
objectIndex{1} = earth('name','Earth'); % Will be initiallised at the origin.

% ////////////////// DEFINE THE MOONS'S ORBITAL PROPERTIES ////////////////
objectIndex{2} = moon('name','Moon');
orbitalInclination = deg2rad(objectIndex{2}.inclination);
axisProjection = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);
% GENERATE THE MOON IN ORBIT
Moon_config = SBinstance.planarRing(...
    'objects',1,...
    'radius',objectIndex{2}.orbit,...
    'pointA',axisProjection,...
    'pointB',zeros(3,1));
objectIndex{2}.VIRTUAL.globalVelocity = [0;objectIndex{2}.orbitalSpeed;0]; % Initialise with tangential oribit speed of 7.67km/s
objectIndex{2}.VIRTUAL.globalPosition = Moon_config.positions;             % Initialise at the apogee (403k), perigee is 406km
objectIndex{2}.VIRTUAL.quaternion = Moon_config.quaternions;

% ////////////////// DEFINE THE ISS'S ORBITAL PROPERTIES //////////////////
% GENERATE THE ISS IN ORBIT
objectIndex{3} = ISS('name','ISS','scale',inputConfig.scale);
orbitalInclination = deg2rad(0);
axisProjection = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);
% ISS scenario
ISS_config = SBinstance.planarRing(...
    'objects',1,...
    'radius',objectIndex{3}.orbit,...
    'pointA',axisProjection,...
    'pointB',zeros(3,1));
objectIndex{3}.VIRTUAL.globalVelocity = [0;objectIndex{3}.orbitalSpeed;0]; % Initialise with tangential oribit speed of 7.67km/s
objectIndex{3}.VIRTUAL.globalPosition = ISS_config.positions;              % Initialise at the apogee (403k), perigee is 406km
objectIndex{3}.VIRTUAL.quaternion = ISS_config.quaternions;
% PLOT THE SCENE
if inputConfig.plot
    SBinstance.plotObjectIndex(objectIndex);
end
% CLEAR THE REMAINING VARIABLES
clearvars -except objectIndex
end