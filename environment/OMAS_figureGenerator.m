%% OPENMAS FIGURE GENERATOR (OMAS_figureGenerator.m) %%%%%%%%%%%%%%%%
% This script contains an index of the figures that can be requested from
% the simulation using the SIM.figures attribute.

% Author: James A. Douthwaite 23/11/2018

%% INDEX FUNCTION
function [figureNumber] = OMAS_figureGenerator(SIM,objectIndex,DATA,figureNumber,figureLabel)
% INPUTS:
% DATA         - The simulation output DATA structure
% figureNumber - The current figure number
% figureLabel  - The figure identifier
% OUTPUT:
% figureNumber - The updated figure number

% PARAMETERISE ALL THE FIGURES BY CURRENT CONFIGURATION
[DATA.figureProperties] = OMAS_figureProperties(SIM,DATA);
% GENERATE THE REDUCED STATE TRAJECTORY FILE
GetTrajectoryTempFile(SIM,DATA);
% DETERMINE WHICH FIGURE IS TO BE GENERATED
switch upper(char(figureLabel))
    case 'ALL'
        fprintf('[%s]\tAll output figures requested.\n',SIM.phase);
        % MOVE THROUGH THE COMPLETE FIGURE VECTOR
        figureVector = {'EVENTS','COLLISIONS','TRAJECTORIES','SEPARATIONS',...
                        'CLOSEST','INPUTS','PLAN','FIG','GIF','AVI','TIMES'};
        for fig = 1:length(figureVector)
            [figureNumber] = OMAS_figureGenerator(SIM,objectIndex,DATA,figureNumber,figureVector{fig});
        end
        close all;  % Kill figures
        fprintf('[%s]\n[%s]\tAll figures pushed to output directory: \n[%s]\t%s',...
                SIM.phase,SIM.phase,SIM.phase,SIM.outputPath);

    case 'EVENTS'
        fprintf('[%s]\tGenerating the event overview figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_EventOverview(SIM,DATA,figureNumber);  
        
    case 'COLLISIONS'
        fprintf('[%s]\tGenerating collision overview.\n',SIM.phase);
        figureSet = [];
        % CHECK COLLISIONS OCCURED
        if ~isfield(DATA,'uniqueCollisions') || isempty(DATA.uniqueCollisions)
            warning('[%s]\t...No collision data available.\n',SIM.phase);
            return
        end
        for collisionNumber = 1:numel(DATA.uniqueCollisions)
            [figureNumber,figureSet(collisionNumber)] = GetFigure_ObjectCollision(SIM,objectIndex,DATA,figureNumber,DATA.uniqueCollisions(collisionNumber));
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Collision Overview');
        set(windowHandle,'Position',DATA.figureProperties.windowSettings);        % Maximise the figure in the tab
        savefig(windowHandle,[SIM.outputPath,'collision_overview']);              % Save the output figure
        
    case 'TRAJECTORIES'
        fprintf('[%s]\tGenerating global trajectory figure.\n',SIM.phase);
        figureSet = [];
        for objectNum = 1:DATA.totalObjects
            [figureNumber,figureSet(objectNum)] = GetFigure_ObjectTrajectory(SIM,DATA,figureNumber,objectNum);
        end
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Trajectory Overview');
        set(windowHandle,'Position',DATA.figureProperties.windowSettings);        % Maximise the figure in the tab
        savefig(windowHandle,[SIM.outputPath,'globalTrajectory_overview']);       % Save the output figure
        
    case 'SEPARATIONS' 
        fprintf('[%s]\tGenerating object separation figure(s).\n',SIM.phase);
        figureSet = [];
        
        %if SIM.totalAgents < 1
        if (SIM.totalObjects - SIM.totalWaypoints) < 2
            warning('There must be at least two collidable objects to plot seperation data.\n');
            return
        end
        
        for objectNum = 1:(SIM.totalObjects - SIM.totalWaypoints)
            [figureNumber,figureSet(objectNum)] = GetFigure_ObjectSeparations(SIM,DATA,figureNumber,objectNum);
        end
        
        % ASSEMBLE TABBED FIGURE
        windowHandle = GetTabbedFigure(figureSet,'OpenMAS Separation Overview');
        set(windowHandle,'Position', DATA.figureProperties.windowSettings);       % Maximise the figure in the tab
        savefig(windowHandle,[SIM.outputPath,'separations_overview']);            % Save the output figure
        
    case 'CLOSEST'
        fprintf('[%s]\tGenerating closest separation figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_MinimumSeparations(SIM,DATA,figureNumber);
        
    case 'INPUTS'
        fprintf('[%s}\tGenerating control input figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_ControlInputs(SIM,objectIndex,DATA,figureNumber);
        
    case 'PLAN'
        fprintf('[%s}\tGenerating top-down 2D(plan) figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_TopDownView(SIM,objectIndex,DATA,figureNumber);
        
    case 'FIG' 
        fprintf('[%s]\tGenerating isometric trajectory figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_IsometricFigure(SIM,objectIndex,DATA,figureNumber);

    case 'GIF'
        fprintf('[%s]\tGenerating trajectory gif file.\n',SIM.phase);
        [figureNumber,~] = GetFigure_IsometricGif(SIM,objectIndex,DATA,figureNumber);
        
    case 'AVI'
        fprintf('[%s]\tGenerating trajectory avi file.\n',SIM.phase);
        [figureNumber,~] = GetFigure_IsometricAvi(SIM,objectIndex,DATA,figureNumber);
        
    case 'TIMES'
        fprintf('[%s]\tGenerating computation timeseries figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_ComputationTimes(SIM,objectIndex,DATA,figureNumber);
        
    case 'NONE'
        fprintf('[%s]\tFigure output supressed.\n',SIM.phase);
        
    otherwise
        warning('[WARNING] Ignoring figure request "%s": Figure unknown.',char(figureLabel));
end
end

% PREPARATION FUNCTIONS 
function [filePath] = GetTrajectoryTempFile(SIM,DATA)
% This function generates/appends to a temp file containing the trajectory 
% states to be plotted in annimations. The number of states are reduced
% according to the desired frame rate and are lower resolution than the
% true number of steps in some cases.
% OUTPUT:
% filePath        - The path to the system file
% objectIDx       - The trajectory data reduced to the desired number of frames.
% frameTimeVector - The reduced time vector associated with the frames.

% CREATE A TEMP FILE FOR THE GENERATION OF TRAJECTORY ELEMENTS
numRep = floor(SIM.TIME.endStep/DATA.figureProperties.stepsPerFrame);          % Handle numsteps not divisable by stepsPerFrame
logicalMatrix = logical(zeros(1,numRep*DATA.figureProperties.stepsPerFrame));  % Frame selection matrix

% FOR EACH OBJECT, GET THE STATE SET
outputStructure = struct();
for objectNo = 1:DATA.totalObjects
    % GET THE COMPLETE STATE SET
    completeStateSet = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,...
                                                  SIM.globalIDvector,...
                                                  SIM.OBJECTS(objectNo).objectID,...
                                                  inf);                    % All valid states for the object
    completeStateSet = completeStateSet(:,1:numRep*DATA.figureProperties.stepsPerFrame);
    
    % GET THE INDICES OF EACH FRAME IN THE STATE SPACE
    for num = 1:numRep
        logicalMatrix(1,num*DATA.figureProperties.stepsPerFrame) = logical(true);
    end
    % SELECT THE STATES FROM THE COMPLETE SET
    statesToPlot = completeStateSet(:,logicalMatrix);
    % GENERATE DATA STRUCTURE TO BE SAVED
    outputStructure.(sprintf('objectID%d',SIM.OBJECTS(objectNo).objectID)) = statesToPlot;
end
% GET THE TIME VECTOR FOR THE FRAME SET
frameTimeVector = SIM.TIME.timeVector(1:numRep*DATA.figureProperties.stepsPerFrame);
outputStructure.frameTimeVector = frameTimeVector(logicalMatrix);

% SAVE A TEMPORARY FILE WITH THE STATES OF EACH OBJECT
save([SIM.outputPath,SIM.systemFile],'-struct','outputStructure','-append');
end