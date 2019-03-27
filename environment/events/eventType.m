%% EVENT ENUMERATION CLASSES (eventType.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the enumeration definition of the available events
% that can happen inside the simulation environment.

% Author: James A. Douthwaite 05/10/16

classdef eventType < uint8
   % Assign uint8 type to enumeration values
   enumeration
      event(0)          % Generic event
      detection(1)      % Detection
      warning(2)        % Proximity warning 
      collision(3)      % Collision 
      waypoint(4)       % Waypoint achieve      
      null_detection(5) % Detection loss 
      null_warning(6)   % Warning-clear 
      null_collision(7) % Collision-clear
      null_waypoint(8)  % Waypoint drop
   end
end