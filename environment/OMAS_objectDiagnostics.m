%% OpenMAS Object Class Diagnotistic tool (OMAS_objectDiagnostics.m) %%%%%%
% This function is designed to preform a batch analysis of the complete
% class library to verify their development status for simulation.

function [issueCount] = OMAS_objectDiagnostics(varargin)

% Input sanity check
assert(numel(varargin) <= 2,'Please provide an object to test, or nothing to test all objects in the "objects" directory.');

% Generate the path to the model directory
phaseStr = 'DIAGNOSTICS';
outputFileName = 'object-diagnostics.xlsx';

repoPath = mfilename('fullpath');           % Get the system paths
ind = strfind(repoPath,'env');
repoPath = repoPath(1:(ind-1)); 

objectsDir = [repoPath,'objects'];
outputFile = [repoPath,outputFileName];
tempDir = OMAS_system.GetOSPathString([repoPath,'temp']);

% DETERMINE THE INPUTS
switch numel(varargin)
    case 1
        % Assume named model
        targetFiles = varargin(1);
    case 0
        % Assume all models
        targetFiles = dir([objectsDir,'\*.m']);
        targetFiles = {targetFiles.name}';
    otherwise
        error('Please provide a model name, or no input.');
end

fprintf('[%s]\t%d object behaviours found, evaluating...\n',phaseStr,numel(targetFiles));

% Handle file I/O
fclose('all');                                                 % Release any file handle
if exist(outputFile,'file')
    delete(outputFile);                              % Delete output file
end

% Prepare the diagnostic header
headerArray = {'model #','file name','class',...
    'build status','build message',...
    'OMAS status','OMAS message',...
    'notes'};

% Check if file is active
try
    % Try to write the header to the file
    xlswrite(outputFile,headerArray);
catch writeError
    warning('Unable to write to diagnostic file, is it open?');
    rethrow(writeError);
end

% Container for model error status
statusArray  = false(numel(targetFiles),1);
summaryArray = cell(numel(targetFiles),numel(headerArray));

% ///////////////// MOVE THROUGH THE MODELS ///////////////////
for i = 1:numel(targetFiles)
    % For clarity
    fprintf('[%s]\tEvaluating model "%s"...\n',phaseStr,targetFiles{i});
    % Preform diagnostic routine on selected model
    [entry,statusArray(i)] = GetObjectEvaluation(phaseStr,tempDir,targetFiles{i});
    % For clarity
    fprintf('[%s]\tEvaluation for model "%s" complete.\n\n',phaseStr,targetFiles{i});
    
    close all; % Terminate all open figures
        
    summaryArray(i,:) = [num2str(i),entry];                    % Append the test number
end
% /////////////////////////////////////////////////////////////
issueCount = sum(statusArray);                                 % Simply confirm if there are any issues

% Write data to the output .xlsx file
xlswrite(outputFile,[headerArray;summaryArray]);
fclose('all');                                                 % Release any file handle

% Delete temporary directory
try
    rmdir(tempDir,'s');
catch
    warning('There was a problem deleting the temporary file.');
end

% Pass statistics
errorNumber   = sum(statusArray == 0);
successNumber = sum(statusArray == 1);

% Summary string
fprintf('[%s]\tDiagnostics complete, %d object errors, %d objects working correctly.\n',phaseStr,errorNumber,successNumber);
end

% Compute the evaluation of the provided object file
function [entryArray,isSuccessful] = GetObjectEvaluation(phaseStr,tempDir,fileName)

% Input sanity check
assert(ischar(fileName),'The model name must be given as a string.');

% Remove the .m from the file-name for evaluation
if contains(fileName,'.m')
    fileName = strrep(fileName,'.m',''); % Remove extension if necessary
end

% Containers
sim_steps   = 5;
sim_dt      = 0.25;
failString = 'FAILED';
passString = 'PASSED';
importStatus = 0;
simStatus = 0;
isSuccessful = 0;

% Create row for output .csv
entryArray = repmat({'-'},1,7);
entryArray{1} = fileName;
entryArray{2} = 'N/a';

% //////////// ATTEMPT TO INSTANTIATE THE MODEL ///////////////
% CREATE THE OBJECT SET
objectNumber = 2;
objectIndex = cell(objectNumber,1);
try
    for j = 1:objectNumber
        objectIndex{j} = eval(fileName);                                % Attempt to construct the object
    end
    importStatus = 1;
catch instantiationError
    warning(instantiationError.message);
    entryArray{3} = failString;
    entryArray{4} = instantiationError.message;
end
% /////////////////////////////////////////////////////////////

% Check instantiation was successful
if ~importStatus
    fprintf('[%s]\t%s(1) - Object "%s" failed to be instantiated.\n',phaseStr,failString,fileName);
    return
end
fprintf('[%s]\t%s(1) - Object "%s" instantiated.\n',phaseStr,passString,fileName);
 
% Lift basic model parameters
entryArray{2} = class(objectIndex{1});
entryArray{3} = passString;

% ////////// IF THE OBJECT CAN BE SIMULATED, ATTEMPT SIMULATION ///////////
fprintf('[%s]\tSimulating object "%s"\n',phaseStr,fileName);
try
    OMAS_initialise(...
        'objects',objectIndex,...
        'duration',sim_steps*sim_dt,...
        'dt',sim_dt,...
        'idleTimeOut',0.5,...
        'figures','None',...
        'verbosity',0,...
        'outputPath',tempDir);
    simStatus = 1;
catch simError
    warning(simError.message);
    entryArray{5} = failString;
    entryArray{6} = simError.message;
    notes = strcat('[function] -> ',simError.stack(1).name);
    entryArray{7} = strcat(notes,' [line] -> ',num2str(simError.stack(1).line));
end
% /////////////////////////////////////////////////////////////////////////
    
% Check simulation was successful
if ~simStatus
	fprintf('[%s]\t%s(2) - Model "%s" failed simulation test\n.',phaseStr,failString,fileName);
    return
end
fprintf('[%s]\t%s(2) - Model "%s" simulated successfully.\n',phaseStr,passString,fileName);

% Indicate simulation was successful
entryArray{5} = passString;

% Indicate that if there where no errors
isSuccessful = 1; 
end