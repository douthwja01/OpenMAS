%% OPENMAS OBJECT ENUMERATION CLASSES (OMAS_objectType.m) %%%%%%%%%%%%%%%%%
% This file contains the enumeration definition of the simulation
% recognised object types.  

% Author: James A. Douthwaite 21/09/17

classdef OMAS_objectType < uint8
   % Assign uint8 type to enumeration values
   enumeration
      misc(0)           % Generic 
      agent(1)          % Active agent
      obstacle(2)       % Obstacle 
      waypoint(3)       % Waypoints 
   end
end