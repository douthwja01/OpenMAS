%% OPENMAS ANALYSIS/POST METAULATION PROCESSING (OMAS_analyis.m) %%%%%%%%%%%
% This function is designed to prepare the collected METAulation data, 
% conduct some preliminary analysis and object statistics before being 
% output for further analysis. 

% Author: James A. Douthwaite 19/12/2017

function [DATA] = OMAS_analysis(META,objectIndex,EVENTS,DATA)
% INPUT:
% META              - The meta data structure
% objectIndex       - The entity vector
% EVENTS            - The complete event history 
% DATA              - The complete output data structure
% .outputpath       - The path to the output location
% .timevector       - The METAulation time vector
% .globalTrajectories - The system timeseries data [(agents*states) by numsteps]

% OUTPUT:
% DATA              - The comprehensive DATA output data structure

global plotnum

%% INPUT HANDLING
% DETERMINE PLOT PROPERTIES
if ~exist('plotnum','var') || isempty(plotnum)
    plotnum = 1;    % Default to first plot
end

fprintf('[%s]\tOMAS ANALYSIS TOOL...\n[%s]\n',META.phase,META.phase);

% CONFIRM INPUT DATA STRUCTURE
if ~exist('DATA','var') || ~isstruct(DATA)
    warning('The output data structure invalid.')
    return
elseif ~isfield(DATA,'timeVector')
    warning('Required time vector not valid in the output structure.');
    DATA = 0;
    return
end

% SAVE RAW DATA TO DIRECTORY BEFORE PROCESSING
sendOutputToFiles(META,EVENTS,DATA,objectIndex)

%% GET EVENT HISTORY STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to event data parser...\n',META.phase);
try
    % HAND THE EVENT SET TO THE EVENT STATISTICS FUNCTION
    [eventStatistics] = OMAS_eventStatistics(META,EVENTS);            % Collect all events into a subheading
    % COPY EVENT DATA TO OUTPUT STRUCTURE
    if isstruct(eventStatistics)                                           % If event data is present 
        for name = fieldnames(eventStatistics)'                            % Move the event history data to the DATA structure
            DATA.(name{1}) = eventStatistics.(name{1});
        end
    end
catch agentParseError
    warning('[ERROR] A problem occurred parsing the event history data.');
    warning(agentParseError.message);
    fprintf('\n');
end

%% GET AGENT-SIDE STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to agent data parser...\n',META.phase);
try
    [objectIndex,DATA.MEANS] = OMAS_agentStatistics(META,objectIndex);
    DATA.objectIndex = objectIndex;
catch agentParseError
    warning('[ERROR] A problem occurred parsing the agent data.');
    warning(agentParseError.message);
    fprintf('\n');
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

% GET AXIS PROPERTIES RELATING TO THE TRAJECTORY DATA
[minimumStateValues,maximumStateValues] = getTrajectoryAxisProperties(DATA);
figureProperties = struct('axisMinimums',1.1*minimumStateValues,...
                          'axisMaximums',1.1*maximumStateValues,...        % Add 10% for padding
                          'minPosition',min(1.1*minimumStateValues(1:3)),...
                          'maxPosition',max(1.1*maximumStateValues(1:3)));
% ASSIGN THE FIGURE PROPERTIES STRUCTURE TO THE OUTPUT DATA
DATA.figureProperties = figureProperties;

% UPDATE RECORDS OF THE DATA STRUCTURES IN THE OUTPUT FILES
sendOutputToFiles(META,EVENTS,DATA,objectIndex)

% ///////////////////// FIGURE GENERATION PROCEDURE ///////////////////////
% GENERATE OUTPUT FIGURES
for figNum = 1:length(figureList)
%     try
        [plotnum] = OMAS_figureGenerator(META,DATA,plotnum,figureList(figNum)); % Jump to the figure index
%     catch figureGenerationError
%         warning('[ERROR] A problem occurred generating the output figures.');
%         warning(figureGenerationError.message);
%         fprintf('\n');
%     end
end
end

%% RAW DATA HANDLING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET AXIS/STATE LIMITS FROM GLOBAL TRAJECTORY DATA
function [stateMinimums,stateMaximums] = getTrajectoryAxisProperties(DATA)
% This function determines the axis limits for the state vector by
% determining the minimum and maximum values of all states across all
% objects.

% INPUTS:
% DATA           - The METAulation DATA structure
% - totalObjects - The number of agents METAulated
% OUTPUTS:
% stateMinimums - Trajectory axis minimum limit vector
% stateMaximums - Trajectory axis maximum limit vector

% CREATE AXIS DATA CONTAINERS
systemStates = size(DATA.globalTrajectories,1)/DATA.totalObjects;          % The number of states per agent
stateMinimums = zeros(systemStates,1);                   
stateMaximums = zeros(systemStates,1);                                     % Define the trajectory axis limit containers

for ID1 = 1:DATA.totalObjects
    % GET THE AGENTS STATE-TIMESERIES DATA
    [objectStateData] = OMAS_getTrajectoryData(DATA,ID1);                  % Get state timeseries data
    
    % DETERMINE GREATEST MIN/MAX STATES FOR EACH AGENTS
    for x = 1:systemStates
       stateMinimum(x) = min(objectStateData(x,:),[],2);                   % Get the minimum state value for this agent
       stateMaximum(x) = max(objectStateData(x,:),[],2);                   % Get the max state value 
       
       % AFFIRM MINIMUM STATE RECORD
       if stateMinimums(x) > stateMinimum(x)
           stateMinimums(x) = stateMinimum(x);                             % Re-affirm the negative limits
       end
       % AFFIRM MAXIMUM STATE RECORD
       if stateMaximums(x) < stateMaximum(x)
           stateMaximums(x) = stateMaximum(x);                             % Re-affirm the positive limits
       end
    end  
end

% CHECK THAT THE LIMITS ARE VALID
for x = 1:systemStates
    if stateMinimums(x) == stateMaximums(x)
        stateMinimums(x) = stateMinimums(x) - 1;
        stateMaximums(x) = stateMaximums(x) + 1;                           % Modify the axis limits if they match
    end
end

clearvars ID1 k axismax axismin
end

% SAVE DATA TO FILES
function sendOutputToFiles(META,EVENTS,DATA,objectIndex)
% This function is designed to handle the output data from the METAulation
% and export the variables to background files
% INPUTS:
% META        - Local copy of the META variable
% EVENTS      - The cell array of event history objects
% DATA        - The output data structure
% objectIndex - The object (non-META) object class cell array

% DISPLAY ALL VARIABLES
% whos;
% SAVE META DATA
save(strcat(META.outputPath,'META.mat'),'META');
save(strcat(META.outputPath,'OBJECTS.mat'),'objectIndex');
% SAVE EVENT HISTORY
save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
% SAVE OUTPUT DATA
save(strcat(META.outputPath,'DATA.mat'),'DATA');
% DISPLAY NOTIFICATION
% fprintf('[%s]\tData objects outputted to file:\n',META.phase);
% fprintf('[%s]\tDirectory: %s\n',META.phase,META.outputPath);
end