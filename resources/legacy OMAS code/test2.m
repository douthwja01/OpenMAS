
close all; clear all;



velocities = 1;
pointA = [1;1;-1]; pointB = [0;0;0];
%pointA = [0;0;-1]; pointB = [0;0;0];
ringradius = 50;
test = scenario(200);
[stateMatrix, test] = test.random();
%[stateMatrix, test] = test.planarRing(pointA,pointB,ringradius,velocities);

% states

test.plot()
