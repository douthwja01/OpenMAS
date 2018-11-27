%% EARTH-OBSTACLE (earth.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Earth class is an iteration of the spheriod obstacle class aimed
% mostly providing a reference for satelite simulation.

% Author: James A. Douthwaite 27/07/2018

classdef earth < obstacle_spheroid
    properties
        radius = 6.371E+06;         % Radius of the earth (m)
        axialRate = 7.2921150E-5;   % Rotational rate about its vertical axes(rad/s)
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = earth(varargin)
            % This function constructs the cuboid obstacle. The object must
            % be imported and represented with a global position and
            % velocity as all other objects are in OMAS.
                        
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@obstacle_spheroid(varargin); 
                        
            % OVERRIDE OBJECT VIRTUAL PROPERTIES
            obj.VIRTUAL.type = OMAS_objectType.obstacle;
            obj.VIRTUAL.symbol = 'o';
            obj.VIRTUAL.radius = obj.radius;                               % The radius of the earth (m)
            obj.VIRTUAL.colour = single([0,0,1]);                          % Must be a row vector of type single
            % IMPORT THE OBJECT'S GEOMETRY IF IT EXISTS
            [obj.GEOMETRY] = OMAS_graphics.scale(obj.GEOMETRY,obj.VIRTUAL.radius);  % Scale

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end 
        
        % ///////////// UPDATE CYCLE FOR A MOVING OBSTACLE ////////////////
        function [obj] = main(obj,TIME,varargin)
            % SIMPLE UPDATE OF LOCAL STATE
            dt = TIME.dt;
            
            % SIMPLY INTEGRATE THE RATES
            [dXdt] = obj.dynamics_simple(obj.localState(1:6,1),obj.localState(7:9,1),obj.localState(10:12,1));
            eulerState = obj.localState;
            eulerState(1:6,1) = obj.localState(1:6,1) + dt*dXdt;
            
            % UPDATE THE GLOBAL PROPERTIES
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,eulerState);
        end
    end
    % /////////////////////// AGENT/OMAS INTERFACES ///////////////////////
    methods
        % INITIALISE THE [x;x_dot]' 3D STATE VECTOR
        function [obj] = initialise_localState(obj,localXYZVelocity,localXYZrotations)
            % The state initialiser must be called 'initialise_localState'
            % and instead calls the 'initialise_3DVelocities' function in
            % this case. 
            [obj] = obj.initialise_3DVelocities(localXYZVelocity,localXYZrotations);
            % ADD A ROTATIONAL RATE
            obj.localState(12) = obj.axialRate; % The earth rotates a constant rate about its z interial axis
        end
    end
end
