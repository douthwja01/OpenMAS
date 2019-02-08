%% OPENMAS EVENT HANDLER (OMAS_eventHandler.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is designed to handle an event occuring between objectA and
% objectB. An EVENT object is created at the current time, based on the two 
% objects and the enumerated event type. 

% Author: James A. Douthwaite 06/10/17

function [EVENT] = OMAS_eventHandler(time_event,metaObjectA,metaObjectB,eventEnumeration)
% INPUTS:
% time_event         - The time the event are occurred
% metaobjectA        - The reference object (entity or meta)
% metaobjectB        - The detected/collision object (entity or meta)
% eventEnumeration   - The type of event indicated by the simulator

% OUTPUTS:
% EVENT              - The generated detection event structure

% INPUT HANDLING
% if ~exist('eventEnumeration','var')
%     error('[ERROR]\tEvent "type" unspecified.\n'); 
% end

% /////////////// GENERATE THE APPROPRIATE EVENT CLASS ////////////////////
% The simulation will generate one of the following events inreponse to the
% request. The result is an event object and the associated updated meta
% object.
switch eventEnumeration
% DETECTION EVENTS    
    case eventType.detection
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Detection notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = detectionEvent(time_event,metaObjectA,metaObjectB,infoString);
% DETECTION-LOSS EVENTS    
    case eventType.null_detection
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Detection-loss notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = detectionEvent(time_event,metaObjectA,metaObjectB,infoString,eventType.null_detection);
 % WARNING/NEAR-MISS EVENTS
    case eventType.warning
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Proximity warning notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = warningEvent(time_event,metaObjectA,metaObjectB,infoString);
% A WARNING/NEAR MISS CONDITION NULLIFICATION
    case eventType.null_warning        
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Proximity warning-clear notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = warningEvent(time_event,metaObjectA,metaObjectB,infoString,eventType.null_warning);
% COLLISION/GEOMETRIC VIOLATION EVENTS
    case eventType.collision
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Collision notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = collisionEvent(time_event,metaObjectA,metaObjectB,infoString);
% COLLISION/GEOMETRIC VIOLATION NULLIFICATION EVENTS        
    case eventType.null_collision
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Collision-clear notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = collisionEvent(time_event,metaObjectA,metaObjectB,infoString,eventType.null_collision);       
% WAYPOINT ACHIEVED EVENTS   
    case eventType.waypoint
        % GENERATE C-FRIENDLY NOTIFICATION STRING
        infoString = char(['Waypoint notification [',metaObjectA.name,':',metaObjectB.name,']']);
        EVENTobj = waypointEvent(time_event,metaObjectA,metaObjectB,infoString);
% WAYPOINT-LOSS EVENTS
    case eventType.null_waypoint
        % A WAYPOINT RESET EVENT (NON-WAYPOINT)
        infoString = sprintf('Waypoint-reset notification [%s:%s]',metaObjectA.name,metaObjectB.name);
        EVENTobj = waypointEvent(time_event,metaObjectA,metaObjectB,infoString,eventType.null_waypoint);
% EVENT UNKNOWN    
    otherwise
        fprintf('[ERROR]\tEvent type not recognised.');
        EVENT = [];
        return
end
% OUTPUT EVENT STRUCTURE
EVENT = []; 
% CONVERT THE EVENT OBJECT INTO A REPORT STRUCTURE
if ~isempty(EVENTobj)
    % GENERATE THE SIMULATION EVENT STRUCTURE
    [EVENT] = event2Struct(EVENTobj);
    % IF VERBOSE, GENERATE NOTIFICATION OF EVENT CREATION
%     if verbosity > 2
%         fprintf('[%s]\tevent created (ID:%d @t=%.2fs):\t%s\n',EVENTobj.name,EVENTobj.eventID,EVENTobj.time,EVENTobj.info);
%     end
end
end

% GENERATE EVENT RECORD
function [EVENT] = event2Struct(EVENTobject)
% This function converts the event object into a structure for the
% simulation log.
% INPUT
% EVENTobject - The EVENT class object
% OUTPUT
% EVENT       - The EVEN log structure

% Append EVENT CLASS TO STRUCTURE FOR REPORTING
EVENT = struct('eventID',EVENTobject.eventID,...
                  'time',EVENTobject.time,...
                  'type',EVENTobject.type,...
                  'name',EVENTobject.name,...
            'objectID_A',EVENTobject.objectID_A,...
                'name_A',EVENTobject.name_A,...
               'state_A',EVENTobject.state_A,...
            'objectID_B',EVENTobject.objectID_B,...
                'name_B',EVENTobject.name_B,...
               'state_B',EVENTobject.state_B,...    
            'separation',EVENTobject.separation,...
                  'info',EVENTobject.info);               
% clear EVENTobject 
end