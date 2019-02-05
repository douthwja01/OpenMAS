%% ADD ADDITIONAL THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('scenarios'); 

clear all; close all;

% IMPORT OBJECT
object = ARdrone();

% GENERAL SETTINGS
figureProperties = struct();
figureProperties.titleString = sprintf('%s[ID-%d]',object.name,object.objectID);
figureProperties.backgroundColor = 'w';
figureProperties.fontWeight = 'bold';
figureProperties.axisFontSize = 12;
figureProperties.titleFontSize = 18;
figureProperties.orientation = [75 44];
figureProperties.axesVisible = 'off';
figureProperties.axesColor = 'none';
figureProperties.edgeAlpha = 0.02;
figureProperties.edgeColour = 'k';
figureProperties.edgeWidth = 0.1;
figureProperties.faceAlpha = 0.05;
figureProperties.faceColour = 'b';
figureProperties.faceLighting = 'gouraud';
% POSITIONING/ALIGNMENT
figureProperties.ROffset = eye(3);
figureProperties.positionOffset = [-0.05;0;0];
figureProperties.R = OMAS_geometry.eulersToRotationMatrix([deg2rad(25);deg2rad(-25);deg2rad(25)]);
figureProperties.triadScale = 1;
figureProperties.triad_ROffset = OMAS_geometry.eulersToRotationMatrix([0;0;0]);
figureProperties.triadPositionOffset = [-1;-1;-1];

%% GENERATE A PLOT OF THE OBJECTS REPRESENTATION
figureHandle = figure('Name',figureProperties.titleString);
% set(figureHandle,'Color','w');                    % Background colour
ax = axes(figureHandle);
hold on;
% Other plot attributes
title(figureProperties.titleString,'fontweight',figureProperties.fontWeight,'fontsize',figureProperties.titleFontSize);
xlabel('x(m)','fontweight',figureProperties.fontWeight,'fontSize',figureProperties.axisFontSize);
ylabel('y(m)','fontweight',figureProperties.fontWeight,'fontSize',figureProperties.axisFontSize);
zlabel('z(m)','fontweight',figureProperties.fontWeight,'fontSize',figureProperties.axisFontSize);
set(ax,'FontSize',figureProperties.axisFontSize,'fontWeight',figureProperties.fontWeight);
set(ax,'Color',figureProperties.axesColor);
axis('tight','equal');
grid on;  box on;
% AXES-ON LOGICAL
set(ax,'visible',figureProperties.axesVisible)
% SET ORIENTATION
view(figureProperties.orientation);
% REPRESENT GEOMETRY AS A PATCH
object.GEOMETRY = OMAS_graphics.normalise(object.GEOMETRY);
objectHandle = patch(ax,...
    'Vertices',(object.GEOMETRY.vertices*figureProperties.ROffset + figureProperties.positionOffset')*figureProperties.R',...
    'Faces',object.GEOMETRY.faces,...
    'FaceColor',figureProperties.faceColour,...
    'EdgeColor',figureProperties.edgeColour,...
    'EdgeAlpha',figureProperties.edgeAlpha,...
    'FaceAlpha',figureProperties.faceAlpha,...
    'FaceLighting',figureProperties.faceLighting,...
    'LineWidth',0.1);

%% ADD ALIGNED TRIAD
triadHandle = OMAS_graphics.drawTriad(figureHandle,zeros(3,1),figureProperties.R*figureProperties.triad_ROffset,figureProperties.triadScale);
set(triadHandle(1,:),'lineWidth',1.2);
set(triadHandle(1,:),'lineStyle','-');
% set(triadHandle(1,:),'Color','k');

%% ADD GLOBAL TRIAD
triadHandle = OMAS_graphics.drawTriad(figureHandle,figureProperties.triadPositionOffset,eye(3),0.3*figureProperties.triadScale);
set(triadHandle(1,:),'lineWidth',1.2);
set(triadHandle(1,:),'lineStyle','-');
% set(triadHandle(1,:),'Color','k');