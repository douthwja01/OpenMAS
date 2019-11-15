%% EARTH'S MOON (moon.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Moon class is an iteration of the spheriod obstacle class aimed
% mostly providing a reference for satelite simulation.

% Author: James A. Douthwaite 09/02/2019

classdef moon < planetoid
    properties

    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function this = moon(varargin)
            % This function constructs the cuboid obstacle. The object must
            % be imported and represented with a global position and
            % velocity as all other objects are in OMAS.
                        
            % Call the super class
            this = this@planetoid(varargin); 
            
            % Assign default 
            this.inclination  = deg2rad(5.145); % Inclination angle
            this.orbit        = 3.844E+08;      % Orbit around the Earth
            this.orbitalSpeed = 1.022E+03;      % Oribital speed around the Earth
            this.mass         = 7.342E+22;      % Mass of the moon (kg)     
            this.axialRate    = 3;              % Rotational rate about its vertical axes (rad/s)
            this.radius       = 1.737E+06;  	% Radius of the moon (m)
            this = this.SetGLOBAL('colour',single([0,0,1]));% Must be a row vector of type single  
            
            % IMPORT THE OBJECT'S GEOMETRY IF IT EXISTS
            [this.GEOMETRY] = OMAS_graphics.scale(this.GEOMETRY,this.radius);  % Scale
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end           
    end
end
