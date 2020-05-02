%% OPENMAS FIGURE GENERATOR (OMAS_figureGenerator.m) %%%%%%%%%%%%%%%%
% This script contains an index of the figures that can be requested from
% the simulation using the SIM.figures attribute.

% Author: James A. Douthwaite 23/11/2018

%% OPENMAS figure generation interface 
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
    case 'ALL' %% /////////////////////////////////////////////////////////
        fprintf('[%s]\tAll output figures requested.\n',SIM.phase);
        % The complete figure set
        figureVector = {...
            'EVENTS',...
            'COLLISIONS',...
            'TRAJECTORIES',...
            'SEPARATIONS',...
            'CLOSEST',...
            'AVOIDANCE',...
            'INPUTS',...
            'PLAN',...
            'ISOMETRIC',...
            'GIF',...
            'AVI',...
            'TIMES'};
        % Move through the figure set
        for fig = 1:length(figureVector)
            [figureNumber] = OMAS_figureGenerator(SIM,objectIndex,DATA,figureNumber,figureVector{fig});
        end
        % Kill figures
        close all;  
        fprintf('[%s]\n[%s]\tAll figures pushed to output directory: \n[%s]\t%s',...
                SIM.phase,SIM.phase,SIM.phase,SIM.outputPath);
    
    case 'EVENTS' %% //////////////////////////////////////////////////////
        fprintf('[%s]\tGenerating the event overview figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_events(SIM,DATA,figureNumber);  
        
    case 'COLLISIONS' %% //////////////////////////////////////////////////
        fprintf('[%s]\tGenerating collision overview figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_objectCollision(SIM,objectIndex,DATA,figureNumber);
        
    case 'TRAJECTORIES' %% ////////////////////////////////////////////////
        fprintf('[%s]\tGenerating global trajectory figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_objectTrajectory(SIM,DATA,figureNumber);

    case 'SEPARATIONS' %% /////////////////////////////////////////////////
        fprintf('[%s]\tGenerating object separation figure(s).\n',SIM.phase);
        [figureNumber,~] = GetFigure_objectSeparations(SIM,DATA,figureNumber);
        
    case 'CLOSEST' %% /////////////////////////////////////////////////////
        fprintf('[%s]\tGenerating closest separation figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_minimumSeparations(SIM,DATA,figureNumber);
        
    case 'AVOIDANCE' %% ///////////////////////////////////////////////////
        fprintf('[%s]\tGenerating object avoidance rates figure(s).\n',SIM.phase);
        [figureNumber,~] = GetFigure_avoidanceRates(SIM,objectIndex,DATA,figureNumber);
        
    case 'INPUTS' %% //////////////////////////////////////////////////////
        fprintf('[%s}\tGenerating control input figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_inputs(SIM,objectIndex,DATA,figureNumber);
        
    case 'PLAN' %% ////////////////////////////////////////////////////////
        fprintf('[%s}\tGenerating top-down 2D(plan) figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_plan(SIM,objectIndex,DATA,figureNumber);
        
    case 'ISOMETRIC' %% ///////////////////////////////////////////////////
        fprintf('[%s]\tGenerating isometric trajectory figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_isometric(SIM,objectIndex,DATA,figureNumber);

    case 'GIF' %% /////////////////////////////////////////////////////////
        fprintf('[%s]\tGenerating trajectory gif file.\n',SIM.phase);
        [figureNumber,~] = GetFigure_isometricGif(SIM,objectIndex,DATA,figureNumber);
        
    case 'AVI' %% /////////////////////////////////////////////////////////
        fprintf('[%s]\tGenerating trajectory avi file.\n',SIM.phase);
        [figureNumber,~] = GetFigure_isometricAvi(SIM,objectIndex,DATA,figureNumber);
        
    case 'TIMES' %% ///////////////////////////////////////////////////////
        fprintf('[%s]\tGenerating computation timeseries figure.\n',SIM.phase);
        [figureNumber,~] = GetFigure_computationTimes(SIM,objectIndex,DATA,figureNumber);
        
    case 'NONE' %% ////////////////////////////////////////////////////////
        fprintf('[%s]\tFigure output supressed.\n',SIM.phase);
        
    otherwise %% //////////////////////////////////////////////////////////
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
numRep = floor(single(SIM.TIME.endStep)/DATA.figureProperties.stepsPerFrame);          % Handle numsteps not divisable by stepsPerFrame
logicalMatrix = false(1,numRep*DATA.figureProperties.stepsPerFrame);  % Frame selection matrix

% FOR EACH OBJECT, GET THE STATE SET
outputStructure = struct();
for objectNo = 1:DATA.totalObjects
    % GET THE COMPLETE STATE SET
%     completeStateSet = OMAS_getTrajectoryData_mex(DATA.globalTrajectories,...
%                                                   SIM.globalIDvector,...
%                                                   SIM.OBJECTS(objectNo).objectID,...
%                                                   inf);                    % All valid states for the object 
    completeStateSet = OMAS_getTrajectoryData(DATA.globalTrajectories,...
                                              SIM.globalIDvector,...
                                              SIM.OBJECTS(objectNo).objectID,...
                                              inf);                    % All valid states for the object 
    
    completeStateSet = completeStateSet(:,1:numRep*DATA.figureProperties.stepsPerFrame);
    
    % GET THE INDICES OF EACH FRAME IN THE STATE SPACE
    for num = 1:numRep
        logicalMatrix(1,num*DATA.figureProperties.stepsPerFrame) = true;
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
filePath = [SIM.outputPath,SIM.systemFile];
save(filePath,'-struct','outputStructure','-append');
end