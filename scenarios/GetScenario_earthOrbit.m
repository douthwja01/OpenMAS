function [ objectIndex ] = GetScenario_earthOrbit(varargin)
% This function is designed to generate a scenario that represents a
% satellite orbiting earth.

fprintf('[SCENARIO]\tGetting the Earth and satillite scenario.\n');

% DEFAULT CONFIGURATION
defaultConfig = struct(...
    'file','scenario.mat',...
    'agents',[],...
    'agentOrbit',435E3,...       % The true orbit of the ISS
    'agentVelocity',0,...
    'scale',5E3,...
    'waypointOrbit',[],...
    'waypointOffsetAngle',[],...
    'waypointRadius',0.5,...
    'noiseFactor',0,...
    'plot',false);  
                   
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
% Set the global state
objectIndex{2}.SetGLOBAL('position',Moon_config.positions);             % Initialise at the apogee (403k), perigee is 406km
objectIndex{2}.SetGLOBAL('velocity',[0;objectIndex{2}.orbitalSpeed;0]); % Initialise with tangential oribit speed of 7.67km/s
objectIndex{2}.SetGLOBAL('quaternion',Moon_config.quaternions);

% ////////////////// DEFINE THE ISS'S ORBITAL PROPERTIES //////////////////
% GENERATE THE ISS IN ORBIT
objectIndex{3} = ISS('name','ISS','scale',inputConfig.scale);
orbitalInclination = deg2rad(0);
axisProjection    = zeros(3,1);
axisProjection(1) = -1*sin(orbitalInclination);
axisProjection(3) = -1*cos(orbitalInclination);
% ISS scenario
ISS_config = SBinstance.planarRing(...
    'objects',1,...
    'radius',objectIndex{3}.orbit,...
    'pointA',axisProjection,...
    'pointB',zeros(3,1));
% Set the global states
objectIndex{3}.SetGLOBAL('velocity',[0;objectIndex{3}.orbitalSpeed;0]); % Initialise with tangential oribit speed of 7.67km/s
objectIndex{3}.SetGLOBAL('position',ISS_config.positions);              % Initialise at the apogee (403k), perigee is 406km
objectIndex{3}.SetGLOBAL('quaternion',ISS_config.quaternions);
% Plot the scene
if inputConfig.plot           
    SBinstance.plotObjectIndex(objectIndex);
end
clearvars -except objectIndex % Clean-up
end