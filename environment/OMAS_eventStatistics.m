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
% GET COLLISION STATISTICS
[EVENTstatistics.collisions,EVENTstatistics.collisionPercentage] = getEventSummary(SIM,parsedEVENTS.collision);
% GET WAYPOINT STATISTICS
[EVENTstatistics.waypointsAchieved,EVENTstatistics.waypointPercentage] = getEventSummary(SIM,parsedEVENTS.waypoint);
% APPEND THE PARSED EVENT STRUCTURE FOR EXTERNAL MANIPULATION
EVENTstatistics.events = parsedEVENTS;

fprintf('[%s]\t... EVENT statistics collected.\n',SIM.phase);
fprintf('[%s]\n',SIM.phase);

% CLEAR UNUSED VARIABLES
clearvars -except EVENTstatistics 
end

% //////////////////////// EVENTS UTILITIES ///////////////////////////////
% GET STATISTICS OF A GIVEN EVENT TYPE
function [occurances,occurancePercentage] = getEventSummary(SIM,specificEventData)
% This function defines some of the additional collision/performance based
% statistics of the simulation based on the events that have occurred.
% INPUTS:
% SIM                 - The META structure
% EVENTS              - The simulation EVENT history
% OUTPUTS:
% occurances          - The number of event occurances
% occurancePercentage - The percentage of all agents for whom this event occurred 

%% PERCENTAGE SUCCESS DEFINITION
% A problem exists that a trigger may not be processed when for one agent
% if the event occurres exactly between timesteps. i.e temporal resolution
% can cause events not to be triggered/missed when a collision is still
% valid. We instead search for any reference to the ID in the collision
% event history.

occurances = 0;
occurancePercentage = 0;

if isempty(specificEventData)
    return
end
fprintf('[%s]\n',SIM.phase); % Seperate the summaries
fprintf('[%s]\tGetting %s event summary:\n',SIM.phase,specificEventData(1).name);

% MOVE THROUGH THE SIMULATION OBJECT SET
for ID = 1:length(SIM.OBJECTS)
    % ONLY INTERESTED IN EVENTS INVOLVING AGENTS CURRENTLY.
    if SIM.OBJECTS(ID).type == OMAS_objectType.agent
        % FIND INSTANCES OF COLLISIONS FOR AGENT
        logicalA = [specificEventData.objectID_A] == SIM.OBJECTS(ID).objectID; % The cause of an event
        logicalB = [specificEventData.objectID_B] == SIM.OBJECTS(ID).objectID; % The recipient of an event      
        % IF EITHER OCCUR
        if any(logicalA) || any(logicalB)
            occurances = occurances + 1;
            occurancePercentage = occurancePercentage + 100*(1/SIM.totalAgents); % Escrew occurances & percentages
            
            % EXTRACT EVENT
            if isempty(logicalA)
                collisionEvents = specificEventData(logicalB);
            else
                collisionEvents = specificEventData(logicalA); 
            end
            % DISPLAY COLLISION OVERVIEW
            for instance = 1:numel(collisionEvents)
                fprintf('[%s]\t%s(t=%s): %s\twith %s\n',SIM.phase,...
                        upper(collisionEvents(instance).name),num2str(collisionEvents(instance).time),...
                        collisionEvents(instance).name_A,collisionEvents(instance).name_B); 
            end
        end
    end
end
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

