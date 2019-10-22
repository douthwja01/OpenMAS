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
        % CONSTRUCTION METHOD
        function obj = obstacle(varargin)
            % This function is to construct the obstacle object using the
            % object defintions held in the 'objectDefinition' base class.
            
            % Call the super class
            obj = obj@objectDefinition(varargin); 
            
            % Allocate obstacle defaults
            obj = obj.SetVIRTUALparameter('type',OMAS_objectType.obstacle);
            obj = obj.SetVIRTUALparameter('hitBoxType',OMAS_hitBoxType.spherical);
            obj = obj.SetVIRTUALparameter('symbol','o'); 
            obj = obj.SetVIRTUALparameter('priorState',obj.localState); 
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end      
    end
end
%%% STATE VECTOR IS DEFINED AS [x;y;z;psi;the;phi;v;u;w;p;q;r] %%%%%%%%%%%%
