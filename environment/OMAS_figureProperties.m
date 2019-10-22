% /////////////////////// GLOBAL FIGURE PROPERTIES ////////////////////////
function [figureProperties] = OMAS_figureProperties(SIM,DATA)
% This function contains a set of predefined axis parameters for figure
% generation. Prespecified values are used to ensure regularity across all
% output figures.
% INPUT:
% SIM               - A local instance of the META object
% intialProperties  - The previous property structure to update
% OUTPUT:
% figureProperties  - A structure of the figure properties

% /////////////////////// GENERAL FIGURE SETUP ////////////////////////////
figureProperties = struct();
figureProperties.cells = 2;                     % Horizontal cells
figureProperties.alignment = 10;                % Horizontal window alignment
figureProperties.margin = 30;                   % Percentage signal margin of plots
figureProperties.spacing = 40;                  % Spacing between figures
figureProperties.tailLength = 4;                % The tail length (in seconds) of comet-like figures
figureProperties.titleFontSize = 28;
figureProperties.axisFontSize = 18;
figureProperties.fontWeight = 'bold';
figureProperties.fontName = 'Helvetica';
figureProperties.interpreter = 'latex';
figureProperties.markerSize = 10;
figureProperties.markerEdgeColor = 'k';
figureProperties.lineWidth = 2;                 % Applies to both the marker and trails
figureProperties.lineStyle = ':';
figureProperties.edgeColor = 'k';
figureProperties.edgeAlpha = 0.2;     
figureProperties.patchLineWidth = 1;            % Patch properties
figureProperties.faceLighting = 'gouraud';
figureProperties.faceAlpha = 0.7;
figureProperties.figureColor = 'w';
figureProperties.axesColor = 'w';               % Universal background colour grey: [0.9 0.9 0.9]    
figureProperties.publish = false;

% /////////////////////// SCREEN CONFIGURATION ////////////////////////////
set(0,'units','pixels')
figureProperties.screensize = get(0,'ScreenSize');
figureProperties.windowSettings = [figureProperties.alignment; 50;
                                   0.6*figureProperties.screensize(3);
                                   0.8*figureProperties.screensize(4)];   % [x y width height]   

% IF ONLY THE BASIC PROPERTIES ARE NEEDED
if nargin == 0
    return
end

figureProperties.publish = SIM.publishMode;     % Generate the pdfs associated with publication

% ERROR CATCHING
if figureProperties.tailLength < SIM.TIME.dt                               % Only needed for long simulations with large time-steps
    figureProperties.tailLength = SIM.TIME.dt;
end

% OBJECT MAXIMAL GEOMETRIES
figureProperties.objectMaximalRadii = max([SIM.OBJECTS.radius]);

% Augment the axis maximum/miniums
[stateMinimums,stateMaximums] = GetTrajectoryAxisProperties(DATA);
figureProperties.axisMinimums = stateMinimums;
figureProperties.axisMaximums = stateMaximums;
figureProperties.minPosition  = min(figureProperties.axisMinimums(1:3));
figureProperties.maxPosition  = max(figureProperties.axisMaximums(1:3));
if abs(figureProperties.minPosition) < abs(figureProperties.maxPosition)
    figureProperties.maxAbsPosition = abs(figureProperties.maxPosition);
else
    figureProperties.maxAbsPosition = abs(figureProperties.minPosition);
end
                
% PREPARE THE DEFAULT LEGEND SET
figureProperties.legendEntries = cell(length(SIM.OBJECTS),1);
for entry = 1:length(SIM.OBJECTS)
    figureProperties.legendEntries{entry} = sprintf('[ID-%s] %s',num2str(SIM.OBJECTS(entry).objectID),SIM.OBJECTS(entry).name);
end

% ///////////////////// ANNIMATION CONFIGURATIONS /////////////////////////
% CONFIGURE VIDEO SETTINGS
fps = 50;                                                                  % Define a maximum frame rate 
if SIM.TIME.endStep > fps*SIM.TIME.endTime
    numFrames = fps*SIM.TIME.endTime; 
else
    numFrames = SIM.TIME.endStep;                                          % The default number of frames
    fps = numFrames/SIM.TIME.endTime;
end

if fps < 1
    fps = 1;
end
figureProperties.fps = fps;

% CALCULATE THE RESULTANT STEPS PER FRAME
stepsPerFrame = floor(single(SIM.TIME.endStep)/single(numFrames));                 % Get the nearest integer steps/frame
if stepsPerFrame == 0                                                      % If the frame frequency is higher than the simulation sample frequency
    stepsPerFrame = 1;                                                     % Default to one frame per sample (highest sample rate)
end
figureProperties.stepsPerFrame = stepsPerFrame;
end

% GET AXIS/STATE LIMITS FROM GLOBAL TRAJECTORY DATA
function [globalMinimums,globalMaximums] = GetTrajectoryAxisProperties(DATA)
% This function determines the axis limits for the state vector by
% determining the minimum and maximum values of all states across all
% objects.

% INPUTS:
% DATA           - The METAulation DATA structure
% - totalObjects - The number of agents METAulated
% OUTPUTS:
% globalMinimums - Trajectory axis minimum limit vector
% globalMaximums - Trajectory axis maximum limit vector

% CREATE AXIS DATA CONTAINERS
systemStates = size(DATA.globalTrajectories,1)/DATA.totalObjects;          % The number of states per agent
tempMatrix   = zeros(systemStates,2*DATA.totalObjects);
tempIDvector = 1:1:DATA.totalObjects;
for ID1 = 1:DATA.totalObjects
    % GET THE AGENTS STATE-TIMESERIES DATA
	% DONT CARE ABOUT objectIDs AT THIS POINT
    [objectStateData] = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,...
                                                   uint16(tempIDvector),...
                                                   uint16(ID1),...
                                                   inf);
    % GET TIME-SERIES MIN-MAX VALUES
    tempMatrix(:,ID1) = min(objectStateData,[],2);                         % Object state minimums
    tempMatrix(:,systemStates+ID1) = max(objectStateData,[],2);            % Object state minimums
end
% COMPUTE THE MIN/MAX'S FOR THE COMPLETE OBJECT SET
globalMinimums = min(tempMatrix,[],2);
globalMaximums = max(tempMatrix,[],2);  

clearvars ID1 objectStateData tempMatrix
end