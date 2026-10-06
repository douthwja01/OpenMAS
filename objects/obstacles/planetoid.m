%% PLANETOID CLASS (planetoid.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Moon class is an iteration of the spheriod obstacle class aimed
% mostly providing a reference for satelite simulation.

% Author: James A. Douthwaite 09/02/2019

classdef planetoid < obstacle_spheroid
    properties
        inclination;    % Inclination angle
        orbit;          % Orbit 
        orbitalSpeed;   % Oribital speed 
        mass;           % Planetoid mass
        axialRate = 1;  % Rotational rate about its vertical axes(rad/s)
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function this = planetoid(varargin)
            % This function constructs the cuboid obstacle. The object must
            % be imported and represented with a global position and
            % velocity as all other objects are in OMAS.
                        
            % Call the super class
            this = this@obstacle_spheroid(varargin); 
                        
            % IMPORT THE OBJECT'S GEOMETRY IF IT EXISTS
            [this.GEOMETRY] = OMAS_graphics.scale(this.GEOMETRY,this.radius);  % Scale
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end           
        % Setup - X = [x;x_dot]' 3D STATE VECTOR
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % The state initialiser must be called 'initialise_localState'
            % and instead calls the 'initialise_3DVelocities' function in
            % this case. 
            [this] = this.setup_3DVelocities(localXYZVelocity,localXYZrotations);
            % ADD A ROTATIONAL RATE
            this.localState(12) = this.axialRate; % The earth rotates a constant rate about its z interial axis
        end
        % Main
        function [this] = main(this,TIME,varargin)
            % SIMPLE UPDATE OF LOCAL STATE
            dt = TIME.dt;
            X = this.localState(1:6,1);
            U = this.localState(7:12,1);
            
            % SIMPLY INTEGRATE THE RATES
            [dXdt] = this.SingleIntegratorDynamics(X,U);
            eulerState = this.localState;
            eulerState(1:6,1) = this.localState(1:6,1) + dt*dXdt;
            
            % UPDATE THE GLOBAL PROPERTIES
            [this] = this.GlobalUpdate_3DVelocities(dt,eulerState);
        end
    end
end
