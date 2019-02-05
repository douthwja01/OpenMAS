%% GLOBAL HAWK FIXED WING AIRCRAFT (globalHawk.m) /////////////////////////

% Author: James A. Douthwaite

classdef globalHawk < fixedWing
    properties
        % Any specific properties of globalhawk
    end
    %  CLASS METHODS
    methods
        % CONSTRUCTOR METHOD
        function obj = globalHawk(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@fixedWing(varargin);                                       % Create the super class 'agent'
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
    end
end