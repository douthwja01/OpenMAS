% ///////////////////////// PREPARATION FUNCTIONS /////////////////////////
function [filePath] = generateTrajectoryTempFile(SIM,DATA)
% This function generates/appends to a temp file containing the trajectory 
% states to be plotted in annimations. The number of states are reduced
% according to the desired frame rate and are lower resolution than the
% true number of steps in some cases.
% OUTPUT:
% filePath        - The path to the system file
% objectIDx       - The trajectory data reduced to the desired number of frames.
% frameTimeVector - The reduced time vector associated with the frames.

% CREATE A TEMP FILE FOR THE GENERATION OF TRAJECTORY ELEMENTS
numRep = floor(SIM.TIME.numSteps/DATA.figureProperties.stepsPerFrame);          % Handle numsteps not divisable by stepsPerFrame
logicalMatrix = logical(zeros(1,numRep*DATA.figureProperties.stepsPerFrame));   % Frame selection matrix

% FOR EACH OBJECT, GET THE STATE SET
outputStructure = struct();
for indexValue = 1:DATA.totalObjects
    % GET THE COMPLETE STATE SET
    [completeStateSet] = OMAS_getTrajectoryData(DATA,indexValue);          % All valid states for the object
    completeStateSet = completeStateSet(:,1:numRep*DATA.figureProperties.stepsPerFrame);
    % GET THE INDICES OF EACH FRAME IN THE STATE SPACE
    for num = 1:numRep
        logicalMatrix(1,num*DATA.figureProperties.stepsPerFrame) = logical(true);
    end
    % SELECT THE STATES FROM THE COMPLETE SET
    statesToPlot = completeStateSet(:,logicalMatrix);
    % GENERATE DATA STRUCTURE TO BE SAVED
    outputStructure.(sprintf('objectID%d',SIM.OBJECTS(indexValue).objectID)) = statesToPlot;
end
% GET THE TIME VECTOR FOR THE FRAME SET
frameTimeVector = SIM.TIME.timeVector(1:numRep*DATA.figureProperties.stepsPerFrame);
outputStructure.frameTimeVector = frameTimeVector(logicalMatrix);

% SAVE A TEMPORARY FILE WITH THE STATES OF EACH OBJECT
filePath = [SIM.outputPath,SIM.systemFile];
save(filePath,'-struct','outputStructure','-append');
end