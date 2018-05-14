%% BASIC WAYPOINT EVENT (waypointEvent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed as a child class to a generic "event", used to 
% define the characteristics of a 'waypoint achieved' event.

% Author: James A. Douthwaite

classdef waypointEvent < interactionEvent
    
%% INITIALISE THE WAYPOINT EVENT SPECIFIC PARAMETERS
    properties
        % Waypoint unique properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR THE WAYPOINT EVENT
        function obj = waypointEvent(time,objectA,objectB,summaryInfo,enum)
            % INPUTS:
            % time          - Detection time (s)
            % objectA       - Class for the reference object
            % objectIDB     - Class for the second object
            % summaryInfo   - Event description string
            % enum          - Provided alternative enum (nullified?)
            % OUTPUTS:
            % obj           - The waypoint event object
                        
            % INITIALISE THE GENERIC EVENT SUPERCLASS
            obj = obj@interactionEvent(time,'WAYPOINT',objectA,objectB,summaryInfo);
 
            % ASSIGN AN EVENT NUMERATION         
%             if ~exist('enum','var') || ~isa(enum,'uint8')  
            if nargin < 5 || ~isa(enum,'uint8')  
                obj.type = eventType.waypoint;
            else
                obj.type = enum;
            end
        end
    end
end