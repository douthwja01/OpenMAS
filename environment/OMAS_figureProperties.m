% /////////////////////// GLOBAL FIGURE PROPERTIES ////////////////////////
function [figureProperties] = OMAS_figureProperties(metaObjectSet,initialProperties)
% This function contains a set of predefined axis parameters for figure
% generation. Prespecified values are used to ensure regularity across all
% output figures.
% INPUT:
% metaObjectSet     - A vector containing the object meta properties
% intialProperties  - The previous property structure to update
% OUTPUT:
% figureProperties  - A structure of the figure properties

% GENERAL FIGURE SETUP
figureProperties = struct('cells',2,...                                    % Horezontal cells
                      'alignment',10,...                                   % Horezontal window alignment
                         'margin',30,...                                   % Percentage signal margin of plots
                        'spacing',40,...                                   % Spacing between figures
                     'tailLength',4,...                                    % The tail length (in seconds) of comet-like figures
                  'titleFontSize',24,...
                   'axisFontSize',18,...
                     'fontWeight','bold',...
                     'MarkerSize',10,...
                'MarkerEdgeColor','k',...
                      'LineWidth',2,...                                    % Applies to both the marker and trails
                      'LineStyle',':',...
                      'EdgeColor','k',...
                      'EdgeAlpha',0.2,...     
                 'PatchLineWidth',1,...                                    % Patch properties
                   'FaceLighting','gouraud',...
                    'figureColor','w',...
                      'axesColor','w');                                    % Universal background colour grey: [0.9 0.9 0.9]    
               
% APPEND SCREEN-SPECFIC DATA
set(0,'units','pixels')
figureProperties.screensize = get(0,'ScreenSize');
figureProperties.windowSettings  = [figureProperties.alignment; 50;
                                    0.6*figureProperties.screensize(3);
                                    0.8*figureProperties.screensize(4)];   % [x y width height]
                  
% PREPARE THE DEFAULT LEGEND SET
figureProperties.legendEntries = cell(length(metaObjectSet),1);
for entry = 1:length(metaObjectSet)
    figureProperties.legendEntries{entry} = sprintf('[ID:%s] %s',num2str(metaObjectSet(entry).objectID),metaObjectSet(entry).name);
end

% IF AN EXTERNAL PROPERTY STRUCTURE IS PROVIDED, UPDATE IT WITH
% AFOREMENTIONED FIELDS
if nargin == 2
    externalFields = fieldnames(initialProperties);
    for i = 1:numel(externalFields)
        figureProperties.(externalFields{i}) = initialProperties.(externalFields{i});
    end
end

end