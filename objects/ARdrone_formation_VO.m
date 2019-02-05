%% ARdrone WITH FORMATION CONTROL AND COLLISION AVOIDANCE %%%%%%%%%%%%%%%%%
% This agent incorperates 
% Author: James A. Douthwaite

classdef ARdrone_formation_VO < ARdrone_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE ARdrone/formation AGENT HYBRID
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = ARdrone_formation_VO(varargin)
            % This function is to construct the ARdrone object using the
            % object defintions held in the 'agent' class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object

            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@ARdrone_formation(varargin);
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);   
        end
        % OVERRIDE THE MAIN CYCLE
        function [obj] = main(obj,ENV,varargin)
            
            % DO SOMETHING
        end
    end
    % 
    methods
        
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]