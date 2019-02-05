%% BASIC COLLISION EVENT (collisionEvent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed as a child class to a generic "event", used to 
% define the characteristics of a collision event.

% Author: James A. Douthwaite

classdef collisionEvent < interactionEvent
    
%% INITIALISE THE COLLISION EVENT SPECIFIC PARAMETERS
    properties
        % Collision unique properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR THE COLLISION EVENT
        function obj = collisionEvent(time,objectA,objectB,summaryInfo,enum)
            % INPUTS:
            % time          - Detection time (s)
            % objectA       - Class for the reference object
            % objectIDB     - Class for the second object
            % summaryInfo   - Event description string
            % enum          - Provided alternative enum (nullified?)
            % OUTPUTS:
            % obj           - The collision object
                        
            % INITIALISE THE GENERIC EVENT SUPERCLASS
            obj = obj@interactionEvent(time,'COLLISION',objectA,objectB,summaryInfo);
 
            % ASSIGN AN EVENT NUMERATION         
%             if ~exist('enum','var') || ~isa(enum,'uint8')  
            if nargin < 5 || ~isa(enum,'uint8') 
                obj.type = eventType.collision;
            else
                obj.type = enum;
            end
        end
    end
end