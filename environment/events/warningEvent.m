%% WARNING EVENT (detectionEvent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed as a child class to a generic "event", used to 
% define the characteristics of an near-miss event.

% Author: James A. Douthwaite 02/10/2016

classdef warningEvent < interactionEvent
    
%% INITIALISE THE COLLISION EVENT SPECIFIC PARAMETERS
    properties
        % Warning unique properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR THE WARNING EVENT
        function obj = warningEvent(time,objectA,objectB,summaryInfo,enum)
            % INPUTS:
            % time          - Warning time (s)
            % objectIDA     - Object ID for the reference object
            % objectIDB     - Object ID for the second object
            % summaryInfo   - Event description string
            % enum          - Provided alternative enum (nullified?)
            % OUTPUTS:
            % obj           - The detection object
                                   
            % INITIALISE THE GENERIC INTERACTIONS SUPERCLASS
            obj = obj@interactionEvent(time,'WARNING',objectA,objectB,summaryInfo);
            
            % ASSIGN AN EVENT NUMERATION         
%             if ~exist('enum','var') || ~isa(enum,'uint8')  
            if nargin < 5 %~isa(enum,'uint8')
                obj.type = eventType.warning;
            else
                obj.type = enum;
            end
        end
    end
end