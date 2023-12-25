%% ARDRONE FORMATION AGENT (ARdrone_formation.m) %%%%%%%%%%%%%%%%%%%%%%%%%%
% This agent incorperates both the formation tracking elements from the 
% "agent_formation" class and the system plant and control dynamics from
% the "ARdrone_LQR" class.

% Author: James A. Douthwaite

classdef ARdrone_formation < ARdrone_LQR & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE ARdrone/formation AGENT HYBRID
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function obj = ARdrone_formation(varargin)
            % This function is to construct the ARdrone object using the
            % object defintions held in the 'agent' class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object

            % Call the super class
            obj = obj@ARdrone_LQR(varargin);
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]