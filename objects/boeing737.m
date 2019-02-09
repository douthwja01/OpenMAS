%% BOEING 737 FIXED WING AIRCRAFT (boeing737.m) /////////////////////////

% Author: James A. Douthwaite

classdef boeing737 < fixedWing
    properties
        % Any specific properties of globalhawk
    end
    %  CLASS METHODS
    methods
        % CONSTRUCTOR METHOD
        function obj = boeing737(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@fixedWing(varargin);                                       % Create the super class 'agent'
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
    end
end