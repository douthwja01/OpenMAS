%% BASIC DETECTION EVENT (detectionEvent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed as a child class to a generic "event", used to 
% define the characteristics of an object detection event.

% Author: James A. Douthwaite 28/05/2016

classdef detectionEvent < interactionEvent
    
%% INITIALISE THE COLLISION EVENT SPECIFIC PARAMETERS
    properties
        % Detection unique properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR DETECTION
        function obj = detectionEvent(time,objectA,objectB,summaryInfo,enum)
            % INPUTS:
            % time          - Detection time (s)
            % objectIDA     - Object ID for the reference object
            % objectIDB     - Object ID for the second object
            % summaryInfo   - Event description string
            % enum          - Provided alternative enum (nullified?)
            % OUTPUTS:
            % obj           - The detection object
            
            % INITIALISE THE GENERIC INTERACTIONS SUPERCLASS
            obj = obj@interactionEvent(time,'DETECTION',objectA,objectB,summaryInfo);
            
            % ASSIGN AN EVENT NUMERATION          
%             if ~exist('enum','var') || ~isa(enum,'uint8')    	
            if nargin < 5 || ~isa(enum,'uint8') 
                obj.type = eventType.detection;
            else
                obj.type = enum;
            end
        end
    end
end