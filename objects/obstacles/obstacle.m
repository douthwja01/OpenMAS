%% GENERIC OBSTACLE CLASS (obstacle.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic obstacle and import its variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite 20/05/2016

classdef obstacle < objectDefinition
%%% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % OBSTACLE VIRTUAL PROPERTIES
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = obstacle(varargin)
            % This function is to construct the obstacle object using the
            % object defintions held in the 'objectDefinition' base class.
            
            % Call the super class
            this = this@objectDefinition(varargin); 
            
            % Allocate obstacle defaults
            this.localState = zeros(6,1);
            this.SetGLOBAL('type',OMAS_objectType.obstacle);
            this.SetGLOBAL('hitBoxType',OMAS_hitBoxType.spherical);
            this.SetGLOBAL('symbol','o'); 
            this.SetGLOBAL('priorState',this.localState); 
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end      
    end
end
%%% STATE VECTOR IS DEFINED AS [x;y;z;psi;the;phi;v;u;w;p;q;r] %%%%%%%%%%%%
