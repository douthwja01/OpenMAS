%% GLOBAL HAWK FIXED WING AIRCRAFT (globalHawk.m) /////////////////////////

% Author: James A. Douthwaite

classdef globalHawk < fixedWing
    properties
        % Any specific properties of globalhawk
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function obj = globalHawk(varargin)
            % Call the super class
            obj@fixedWing(varargin);   
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
    end
end