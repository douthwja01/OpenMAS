%% OPENMAS EVENT STATISTICS CALCULATOR (OMAS_eventStatics.m) %%%%%%%%%%%%%%
% This function is designed to parse the output event history data
% structure list in order to calculate session event statistics.

% Author: James A. Douthwaite 07/10/2016

function [EVENTstatistics] = OMAS_eventStatistics(SIM,EVENTS)
% Each event cant be identified by its .type and .name
% INPUTS:
% SIM       - The simultaion META structure
% EVENTS    - The simulation event history structure
% OUTPUTS:
% EVENTstatistics - The event statistical data

fprintf('[%s]\tParsing simulation EVENT data...\n',SIM.phase);

% INPUT HANDLING
if isempty(EVENTS)
    fprintf('[%s]\tNo EVENT history to parse.\n',SIM.phase);
    EVENTstatistics = [];
    return
elseif ~isstruct(SIM) || ~isstruct(EVENTS)
    warning('[ERROR] Parse error: the event history is an invalid structure.');
    EVENTstatistics = [];
    return
end

% //////////////////////// GENERAL EVENT DATA /////////////////////////////
% PARSE THE EVENT TYPES
[parsedEVENTS] = parseEventTypes(EVENTS);
% DERIVE THE TIMESERIES DATA
[EVENTstatistics.eventTimeSeries] = getEventTimeSeries(SIM.TIME,parsedEVENTS);
% GET THE TOTAL NUMBER OF EVENTS
EVENTstatistics.totalEvents = numel(EVENTS);                               % The length of the EVENTS output structure

% //////////////////// SPECIFIC EVENT DATA ////////////////////////////////
% >>>>>>> GET COLLISION STATISTICS
[~,uniqueCollisions] = getEventSummary(SIM,parsedEVENTS.collision);
% WE ARE INTERESTED IN THE COLLISIONS RELATING ONLY TO THE AGENTS
if isempty(uniqueCollisions)
    collisionIDs = [];
else
    collisionIDs = unique([uniqueCollisions(:).objectID_A,uniqueCollisions(:).objectID_B]); % The unique IDs involved in collisions
end
EVENTstatistics.uniqueCollisions = uniqueCollisions;

agentSet = SIM.OBJECTS([SIM.OBJECTS.type] == OMAS_objectType.agent);   % The complete agent set
% THE AGENTS THAT COLLIDED/WHERE COLLIDED WITH
collidedAgents = SIM.OBJECTS(ismember([agentSet.objectID],collisionIDs));
% SUMMARISE COLLISIONS
EVENTstatistics.collisions = numel(collidedAgents);
EVENTstatistics.collisionPercentage = 100*(EVENTstatistics.collisions/double(SIM.totalAgents));

% >>>>>>> GET WAYPOINT STATISTICS
[EVENTstatistics.waypointsAchieved,~] = getEventSummary(SIM,parsedEVENTS.waypoint);
EVENTstatistics.waypointPercentage = 100*(EVENTstatistics.waypointsAchieved/double(SIM.totalWaypoints));
if isnan(EVENTstatistics.waypointPercentage)
    EVENTstatistics.waypointPercentage = 0;
end

% APPEND THE PARSED EVENT STRUCTURE FOR EXTERNAL MANIPULATION
EVENTstatistics.events = parsedEVENTS;
% PRINT THE SUMMARY
fprintf('[%s]\n',SIM.phase);
fprintf('[%s]\tEVENT statistics collected.\n',SIM.phase);
fprintf('[%s]\n',SIM.phase);

% CLEAR UNUSED VARIABLES
clearvars -except EVENTstatistics 
end

% //////////////////////// EVENTS UTILITIES ///////////////////////////////
% GET STATISTICS OF A GIVEN EVENT TYPE
function [uniqueOccurances,uniqueEvents] = getEventSummary(SIM,specificEventData)
% This function defines some of the additional collision/performance based
% statistics of the simulation based on the events that have occurred.
% INPUTS:
% SIM                 - The META structure
% EVENTS              - The simulation EVENT history
% OUTPUTS:
% uniqueOccurances    - The number of unique event occurances
% uniqueEvents        - The unique event objects  

% A problem exists that a trigger may not be processed when for one agent
% if the event occurres exactly between timesteps. i.e temporal resolution
% can cause events not to be triggered/missed when a collision is still
% valid. We instead search for any reference to the ID in the collision
% event history.

uniqueOccurances = 0;
uniqueEvents = [];

% INPUT HANDLING
if isempty(specificEventData)
    return
end
assert(numel(specificEventData) == sum([specificEventData.type] == specificEventData(1).type),...
       'Provided event data is not all of the same type.');
   
% BEGIN PARSING EVENT STRUCTURES   
fprintf('[%s]\n',SIM.phase); % Seperate the summaries
fprintf('[%s]\tGetting %s event summary:\n',SIM.phase,specificEventData(1).name);

% Situations may occur where an agent will collide and then continue
% before repeating the collision. We wish to represent this as two
% seperate collision events.

% CONTAINERS
illegalEvents = [];
% PARSE EVENTS FOR UNIQUE OCCURANCES OF THE EVENT/INTERACTION TYPE
for eventNum = 1:numel(specificEventData)
    evaluationEvent = specificEventData(eventNum);
    % GET THE EVENTS WITH THIS AGENT ID AS THE RECIPIENT
    reflectedEvents = specificEventData([specificEventData.objectID_B] == evaluationEvent.objectID_A);
    % GET THE EVENTS WITH THE REFLECTED PROPERTIES
    reciprocalEvents = reflectedEvents([reflectedEvents.objectID_A] == evaluationEvent.objectID_B);
    % DETERMINE IF THE INTERACTION IS ONE WAY 
    if numel(reciprocalEvents) == 0 
        % There is no reciprication, by definition is unique
        uniqueEvents = [uniqueEvents,evaluationEvent];
    else
        for t = 1:numel(reciprocalEvents)
            % OCCURS AT THE SAME TIME
            temporalCondition = evaluationEvent.time == reciprocalEvents(t).time;
            % IS NOT THE REFLECTION OF A PREVIOUS EVENT
            illegalEventCondition = any(ismember(illegalEvents,reciprocalEvents(t).eventID));
            % DETERMINE IF UNIQUE
            if temporalCondition && ~illegalEventCondition
                uniqueEvents = [uniqueEvents,evaluationEvent];
                illegalEvents = [illegalEvents;evaluationEvent.eventID];
            end
        end
    end
end

% DISPLAY OVERVIEW
for instance = 1:numel(uniqueEvents)
    fprintf('[%s]\t%s(t=%.2fs): %s\twith %s\n',SIM.phase,...
        upper(uniqueEvents(instance).name),...
        uniqueEvents(instance).time,...
        uniqueEvents(instance).name_A,...
        uniqueEvents(instance).name_B);
end
% CALCULATE PRELIMINARY STATS
uniqueOccurances = numel(uniqueEvents);
% SHOW THE EVENT SUMMARY
fprintf('[%s]\t ... %s unique %s occurances.\n',SIM.phase,num2str(uniqueOccurances),specificEventData(1).name);
end
% GET THE EVENT TIMESERIES DATA
function [timeSeriesStructure] = getEventTimeSeries(TIME,parsedEVENTS)
% INPUTS:
% SIM             - The simualation overview data
% EVENTdata       - The parsed EVENTS structure
% OUTPUT:
% timeSeriesEvent - The time series data

% Get the event timeseries data
eventStrings = fieldnames(parsedEVENTS);            

% MOVE THROUGH THE EVENT TYPES
for eventName = 1:length(eventStrings)
    eventTimeSeries = zeros(size(TIME.timeVector));                        % Create the empty data containers
    eventStructure = parsedEVENTS.(char(eventStrings(eventName)));         % Get the event substructure
    if numel(eventStructure) ~= 0                                          % No instances, no time-series
        % Convert the number of entries at an indicated times
        for eventNum = 1:length(eventStructure)
            % Get the time the event occurred
            time = eventStructure(eventNum).time;
            indexValue = find(TIME.timeVector == time);                    % Find the index of the time value
            eventTimeSeries(indexValue) = eventTimeSeries(indexValue) + 1;
        end
    end
    % RECORD TIMESERIES FOR EVENT TYPE
    timeSeriesStructure.(char(eventStrings(eventName))) = eventTimeSeries; % Create the new field name
end
clearvars -except timeSeriesStructure
end
% PARSE THE OUTPUT 'EVENT' DATASET
function [summaryData] = parseEventTypes(EVENTS)
% This function returns a structure with the events parsed by their event
% types.

% INPUTS:
% EVENTS    - The structure list of events
% OUTPUTS:
% EVENTdata - The structure with ammended statistics

% GET THE SIMULATION EVENT LISTING
[eventEnum,enumString] = enumeration('eventType');     
for event = 1:numel(eventEnum)
	summaryData.(enumString{event}) = EVENTS([EVENTS.type] == eventEnum(event));
end
end

