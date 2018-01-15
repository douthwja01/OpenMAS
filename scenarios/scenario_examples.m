close all; clear all;

addpath('scenarios');
addpath('environment');

agentNumber = 6;
example_scenario = scenarioBuilder(agentNumber);

%% DEFINE A RING OF 5 OBJECT STATES 
pointA = [0;0;-1]; 
pointB = [0;0.5;0]; 
ringradius = 10; 
velocities = 2;
[config] = example_scenario.planarRing('pointA',pointA,'pointB',pointB,...
                                       'radius',ringradius,'velocities',velocities);
config.plot();

pause;
close all;

%% DEFINE AN ARC OF OBJECT STATES
pointA = [0;0;-1];  
pointB = [0;0;0]; 
ringradius = 10; 
velocities = 4;
offsetAngle = 0.25; % Offset angle in radians
[config] = example_scenario.planarAngle('pointA',pointA,'pointB',pointB,...
                                       'radius',ringradius,'velocities',velocities,...
                                       'offsetAngle',offsetAngle);
config.plot();

pause;
close all;

%% DEFINE A SET OF UNIFORMLY DISTRIBUTED STATES
[config] = example_scenario.randomUniform()
config.plot();

pause;
close all;

%% DEFINE A SET OF NORMALLY DISTRIBUTED STATES
[config] = example_scenario.randomNormal()
% [config] = example_scenario.random()

config.plot();
