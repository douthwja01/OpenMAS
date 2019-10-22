%% INTERACTION EVENT (interactionEvent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This event is designed to be a modifed 'event' to create events that
% resemble interactions between objects and the environment

% Author: James A. Douthwaite 07/10/2016

classdef interactionEvent < eventDefinition

%% INITIALISE THE COLLISION EVENT SPECIFIC PARAMETERS
    properties
        objectID_A = uint32(0);
        name_A;
        state_A;
        objectID_B = uint32(0);
        name_B;
        state_B;
        separation;
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR BASIC INTERACTION EVENTS
        function obj = interactionEvent(time,name,objectA,objectB,summaryInfo)
            % This function initialises the event using the object ID's of both objects involved.
            % INPUTS:
            % time - Detection time (s)
            % objectIDA - Object ID for the reference object
            % objectIDB - Object ID for the second object
            % OUTPUTS:
            % obj - The near-miss object
              
            % INITIALISE THE GENERIC EVENT SUPERCLASS
            obj = obj@eventDefinition(time,name,summaryInfo);       
            
            % DECLARE THE RELEVANT OBJECT DATA
            obj.objectID_A = objectA.objectID;
            obj.objectID_B = objectB.objectID;                             % Object ID's of those involved
            obj.name_A = objectA.name;
            obj.name_B = objectB.name;                                     % Agent names of those involved
            obj.state_A = objectA.globalState;
            obj.state_B = objectB.globalState;                             % Record agent states 
            obj.separation = sqrt(sum((obj.state_B(1:3,1) - obj.state_A(1:3,1)).^2)); % Agent seperation at time of incident
        end
    end
end