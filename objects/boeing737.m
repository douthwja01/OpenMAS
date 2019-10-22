%% BOEING 737 FIXED WING AIRCRAFT (boeing737.m) /////////////////////////

% Author: James A. Douthwaite

classdef boeing737 < fixedWing
    properties
        % Any specific properties of globalhawk
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function obj = boeing737(varargin)
            % Call the super class
            obj@fixedWing(varargin);                                       % Create the super class 'agent'
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
    end
end