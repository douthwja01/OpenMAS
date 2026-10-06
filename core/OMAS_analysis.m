%% OPENMAS ANALYSIS/POST METAULATION PROCESSING (OMAS_analyis.m) %%%%%%%%%%%
% This function is designed to prepare the collected METAulation data, 
% conduct some preliminary analysis and object statistics before being 
% output for further analysis. 

% Author: James A. Douthwaite 19/12/2017

function [DATA] = OMAS_analysis(META,OBJECTS,EVENTS,DATA)
% INPUT:
% META              - The meta data structure
% OBJECTS           - The entity handle vector
% EVENTS            - The complete event history 
% DATA              - The complete output data structure
% .outputpath       - The path to the output location
% .timevector       - The METAulation time vector
% .globalTrajectories - The system timeseries data [(agents*states) by numsteps]

% OUTPUT:
% DATA              - The comprehensive DATA output data structure

% Input sanity check #1
assert(nargin == 4,'Expecting four output structures.');
assert(iscell(OBJECTS),'Expecting a cell array of objects.');

% assert(isstruct(META),'Expecting an OpenMAS META structure');
% assert(isstruct(EVENTS),'Expecting an OpenMAS EVENTS structure');
% assert(isstruct(DATA),'Expecting an OpenMAS DATA structure');

% DETERMINE PLOT PROPERTIES
global plotnum
if ~exist('plotnum','var') || isempty(plotnum)
    plotnum = 1;    % Default to first plot
end

fprintf('[%s]\tOMAS ANALYSIS TOOL...\n[%s]\n',META.phase,META.phase);

% CONFIRM INPUT DATA STRUCTURE
if ~isstruct(DATA)
    warning('The output data structure invalid.')
    return
elseif ~isfield(DATA,'timeVector')
    warning('Required time vector not valid in the output structure.');
    DATA = 0;
    return
end

% SAVE RAW DATA TO DIRECTORY BEFORE PROCESSING
sendOutputToFiles(META,EVENTS,DATA,OBJECTS)

%% GET EVENT HISTORY STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to event data parser...\n',META.phase);
try
    % HAND THE EVENT SET TO THE EVENT STATISTICS FUNCTION
    [eventStatistics] = OMAS_eventStatistics(META,EVENTS);                 % Collect all events into a subheading
    % COPY EVENT DATA TO OUTPUT STRUCTURE
    if isstruct(eventStatistics)                                           % If event data is present 
        for name = fieldnames(eventStatistics)'                            % Move the event history data to the DATA structure
            DATA.(name{1}) = eventStatistics.(name{1});
        end
    end
catch agentParseError
    warning('[OMAS-error] A problem occurred parsing the event history data.');
    warning(agentParseError.message);
    fprintf('\n');
    rethrow(agentParseError);
end

%% GET AGENT-SIDE STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to agent data parser...\n',META.phase);
try
    [OBJECTS,DATA.MEANS] = OMAS_agentStatistics(META,OBJECTS);
catch agentParseError
    warning('[OMAS-error] A problem occurred parsing the agent data.');
    warning(agentParseError.message);
    fprintf('\n');
    rethrow(agentParseError);
end

% //////////////// GENERATE THE DEFAULT FIGURE PARAMETERS /////////////////
fprintf('\n[%s]\tOMAS FIGURE GENERATOR (OMAS_figureGenerator.m).\n[%s]\n',META.phase,META.phase);
% IMPORT THE REQUESTED FIGURE SET FROM THE SESSION DIRECTORY
load([META.outputPath,META.systemFile]); % Import the requested list

% CHECK IF FIGURES HAVE BEEN REQUESTED
if isempty(figureList)
    fprintf('[%s]\tNo output figures requested.\n',META.phase);             % No figures requested by user
    return
end

% UPDATE RECORDS OF THE DATA STRUCTURES IN THE OUTPUT FILES
sendOutputToFiles(META,EVENTS,DATA,OBJECTS)

% ///////////////////// FIGURE GENERATION PROCEDURE ///////////////////////
for figNum = 1:length(figureList)
    try
        % Generate output figures
        [plotnum] = OMAS_figureGenerator(META,OBJECTS,DATA,plotnum,figureList(figNum)); % Jump to the figure index
    catch figureGenerationError
        warning('[OMAS-error] A problem occurred generating the output figures.');
        warning(figureGenerationError.message);
        fprintf('\n');
        rethrow(figureGenerationError);
    end
end
end

% SAVE DATA TO FILES
function sendOutputToFiles(META,EVENTS,DATA,OBJECTS)
% This function is designed to handle the output data from the METAulation
% and export the variables to background files
% INPUTS:
% META        - Local copy of the META variable
% EVENTS      - The cell array of event history objects
% DATA        - The output data structure
% objectIndex - The object (non-META) object class cell array

% DISPLAY ALL VARIABLES
% whos;
save(strcat(META.outputPath,'META.mat'),'META');
save(strcat(META.outputPath,'OBJECTS.mat'),'OBJECTS');
save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
save(strcat(META.outputPath,'DATA.mat'),'DATA');

% DISPLAY NOTIFICATION
% fprintf('[%s]\tData objects outputted to file:\n',META.phase);
% fprintf('[%s]\tDirectory: %s\n',META.phase,META.outputPath);
clearvars; 
end