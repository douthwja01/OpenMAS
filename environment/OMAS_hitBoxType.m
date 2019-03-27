%% OPENMAS HIT-BOX ENUMERATION CLASSES (OMAS_hitBoxType.m) %%%%%%%%%%%%%%%%
% This file contains the enumeration definition of the simulation
% recognised object types.  

% Author: James A. Douthwaite 25/01/19

classdef OMAS_hitBoxType < uint8
   % Assign uint8 type to enumeration values
   enumeration
      none(0);           % Generic 
      spherical(1);      % Spherical bounding box assumption
      AABB(2);           % Axis aligned bounding box
      OBB(3);            % Object orientated bounding box
      capsule(4);        % Capsule enclosed volume
      mesh(5);           % Mesh enclosed volume
   end
end