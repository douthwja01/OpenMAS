%% A10 WARTHOG FIXED WING AIRCRAFT (A10.m) ///////////////////////////////

% Author: James A. Douthwaite

classdef A10 < fixedWing
    properties
        % Any specific properties of globalhawk
    end
    %  CLASS METHODS
    methods
        % CONSTRUCTOR METHOD
        function obj = A10(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@fixedWing(varargin);  % Create the super class 'fixedWing' 
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
    end
end