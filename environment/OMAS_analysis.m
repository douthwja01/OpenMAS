%% OPENMAS ANALYSIS/POST SIMULATION PROCESSING (OMAS_analyis.m) %%%%%%%%%%%
% This function is designed to prepare the collected simulation data, 
% conduct some preliminary analysis and object statistics before being 
% output for further analysis. 

% Author: James A. Douthwaite 19/12/2017

function [DATA] = OMAS_analysis(SIM,objectIndex,EVENTS,DATA)
% INPUT:
% SIM               - The meta data structure
% objectIndex       - The entity vector
% EVENTS            - The complete event history 
% DATA              - The complete output data structure
% .outputpath       - The path to the output location
% .timevector       - The simulation time vector
% .globalTrajectories - The system timeseries data [(agents*states) by numsteps]

% OUTPUT:
% DATA              - The comprehensive DATA output data structure

global plotnum

%% INPUT HANDLING
% DETERMINE PLOT PROPERTIES
if ~exist('plotnum','var') || isempty(plotnum)
    plotnum = 1;    % Default to first plot
end

fprintf('[%s]\tAnalysing object/event data...\n',SIM.phase);

% CONFIRM INPUT DATA STRUCTURE
if ~exist('DATA','var') || ~isstruct(DATA)
    warning('The output data structure invalid.')
    return
elseif ~isfield(DATA,'timeVector')
    warning('Required time vector not valid in the output structure.');
    DATA = 0;
    return
end

%% GET EVENT HISTORY STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to event data parser...\n',SIM.phase);
try
    % HAND THE EVENT SET TO THE EVENT STATISTICS FUNCTION
    [eventStatistics] = OMAS_eventStatistics(SIM,EVENTS);            % Collect all events into a subheading
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
fprintf('[%s]\tMoving to agent data parser...\n',SIM.phase);
try
    [objectIndex,DATA.MEANS] = OMAS_agentStatistics(SIM,objectIndex);
    DATA.objectIndex = objectIndex;
catch agentParseError
    warning('[ERROR] A problem occurred parsing the agent data.');
    warning(agentParseError.message);
    fprintf('\n');
end

%% GENERATE OUTPUT FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\tMoving to figure generator...\n',SIM.phase);
% CHECK IF FIGURES HAVE BEEN REQUESTED
try
    deselectFlag = isempty(SIM.figures) || ~any(SIM.figures); 
    if deselectFlag || strcmpi(SIM.figures,'NONE') 
        fprintf('[%s]\tNo output figures requested.\n',SIM.phase);         % No figures requested by user
        return
    end
catch
end

%% FIGURE GENERATION PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[%s]\tGENERATING OUTPUT FIGURES(OMAS_figureIndex.m).\n',SIM.phase);
% GET GENERIC PREDEFINED AXIS PROPERTIES
[figureProperties] = getGenericAxisProperties(SIM.OBJECTS);                % Get the prespecified axis properties
% GET AXIS PROPERTIES RELATING TO THE TRAJECTORY DATA
[minimumStateValues,maximumStateValues] = getTrajectoryAxisProperties(DATA);
figureProperties.axisMinimums = minimumStateValues;
figureProperties.axisMaximums = maximumStateValues;
figureProperties.minPosition = min(figureProperties.axisMinimums(1:3));
figureProperties.maxPosition = max(figureProperties.axisMaximums(1:3));

% DEFINE THE LEGEND ENTRIES FROM THE EVENT TYPS
[~,enumString] = enumeration('eventType');
for label = 1:numel(enumString)
    figureProperties.legendSet_events{label} = strrep(char(enumString{label}),'_',' ');           % Prepare the label strings from the enumerations
end

% ASSIGN THE FIGURE PROPERTIES STRUCTURE TO THE OUTPUT DATA
DATA.figureProperties = figureProperties;

% GENERATE OUTPUT FIGURES
for figNum = 1:length(SIM.figures)
%     try
        [plotnum] = OMAS_figureGenerator(SIM,DATA,plotnum,SIM.figures(figNum));                              % Jump to the figure index
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
% DATA           - The simulation DATA structure
% - totalObjects - The number of agents simulated
% OUTPUTS:
% stateMinimums - Trajectory axis minimum limit vector
% stateMaximums - Trajectory axis maximum limit vector

% CREATE AXIS DATA CONTAINERS
systemStates = size(DATA.globalTrajectories,1)/DATA.totalObjects;          % The number of states per agent
stateMinimums = zeros(systemStates,1);                   
stateMaximums = zeros(systemStates,1);                                % Define the trajectory axis limit containers

for ID1 = 1:DATA.totalObjects
    % GET THE AGENTS STATE-TIMESERIES DATA
    [objectStateData] = OMAS_getTrajectoryData(DATA,ID1);            % Get state timeseries data
    
    % DETERMINE GREATEST MIN/MAX STATES FOR EACH AGENTS
    for x = 1:systemStates
       stateMinimum(x) = min(objectStateData(x,:),[],2);                       % Get the minimum state value for this agent
       stateMaximum(x) = max(objectStateData(x,:),[],2);                       % Get the max state value 
       
       % AFFIRM MINIMUM STATE RECORD
       if stateMinimums(x) > stateMinimum(x)
           stateMinimums(x) = stateMinimum(x);                            % Re-affirm the negative limits
       end
       % AFFIRM MAXIMUM STATE RECORD
       if stateMaximums(x) < stateMaximum(x)
           stateMaximums(x) = stateMaximum(x);                            % Re-affirm the positive limits
       end
    end  
end

% CHECK THAT THE LIMITS ARE VALID
for x = 1:systemStates
    if stateMinimums(x) == stateMaximums(x)
        stateMinimums(x) = stateMinimums(x) - 1;
        stateMaximums(x) = stateMaximums(x) + 1;                 % Modify the axis limits if they match
    end
end
clearvars ID1 k axismax axismin
end
% OUTPUT FIGURE (GENERIC) PROPERTIES
function [figureProperties] = getGenericAxisProperties(metaObjectSet)
% This function contains a set of predefined axis parameters for figure
% generation. Prespecified values are used to ensure regularity across all
% output figures.
% INPUT:
% metaObjectSet     - A vector containing the object meta properties
% OUTPUT:
% figureProperties  - An structure containing the figure parameters.

%   GENERAL FIGURE SETUP
figureProperties = struct('cells',2,...     % Horezontal cells
                      'alignment',10,...    % Horezontal window alignment
                         'margin',30,...    % Percentage signal margin of plots
                        'spacing',40,...    % Spacing between figures
                     'tailLength',2,...     % The tail length (in seconds) of comet-like figures
                  'titleFontSize',18,...
                   'axisFontSize',10,...
                     'fontWeight','bold',...
                     'MarkerSize',10,...
                'MarkerEdgeColor','k',...
                      'LineWidth',1.5,...     % Applies to both the marker and trails
                      'LineStyle',':');                  
% APPEND SCREEN-SPECFIC DATA
set(0,'units','pixels')
figureProperties.screensize = get(0,'ScreenSize');
figureProperties.size = figureProperties.screensize(1,4)*(4/5);
                  
% PREPARE THE DEFAULT LEGEND SET
figureProperties.legendEntries = cell(length(metaObjectSet),1);
for entry = 1:length(metaObjectSet)
    figureProperties.legendEntries{entry} = sprintf('%s[ID:%s]',metaObjectSet(entry).name,num2str(metaObjectSet(entry).objectID));
end

end