%% EARTH (earth.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Earth class is an iteration of the spheriod obstacle class aimed
% mostly providing a reference for satelite simulation.

% Author: James A. Douthwaite 27/07/2018

classdef earth < planetoid
    properties
        
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function this = earth(varargin)
            % This function constructs the cuboid obstacle. The object must
            % be imported and represented with a global position and
            % velocity as all other objects are in OMAS.
                        
            % CALL THE SUPERCLASS CONSTRUCTOR
            this = this@planetoid(varargin); 

            % Assign "Earth" properties
            this.axialRate = 7.2921150E-5;   % Rotational rate about its vertical axes(rad/s)
            this.radius = 6.371E+06;         % Radius of the earth (m)
            this = this.SetGLOBAL('colour',single([0,0,1]));
            this.GEOMETRY = OMAS_graphics.scale(this.GEOMETRY,this.radius);  % Scale

            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end 
    end
end
