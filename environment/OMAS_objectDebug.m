
function [isIssue] = OMAS_objectDebug(varargin)

% Input sanity check
if numel(varargin) > 1
    error('Please provide an object to test, or nothing to test all objects in the "objects" directory.');
end

% CONSTANTS
[~, userdir] = system('echo %USERPROFILE%'); % Get desktop path
TEST.tempDir = strcat(userdir,'\desktop\temp');
TEST.objectsPath  = '\objects';
TEST.objectNumber = 2;
TEST.phase        = 'DEBUG';
TEST.summaryFile  = [pwd,'\objectDebug.csv'];
isIssue = 0;

% DETERMINE THE INPUTS 
if nargin < 1
    target = dir([cd,TEST.objectsPath,'\*.m']);                % GET ALL .m FILES IN THE OBJECTS DIRECTORIES 
    target = {target.name}'; 
elseif ischar(varargin{1})
    target = varargin(1);
else
    error('Provide the name of the test object.'); 
end

% Write the header to the output file
if exist(TEST.summaryFile,'file')
    delete(TEST.summaryFile);               % Delete 
    delete([strrep(TEST.summaryFile,'.csv',''),'.xlsx']);
end

% GENERATE THE HEADER 
TEST.fid = fopen(TEST.summaryFile,'w');     % Write the string to a file
TEST = getXlsHeader(TEST);
fclose(TEST.fid);                           % Close the summary

% /////////////////////////////////////////////////////////////////////////
% Move through the object
TEST.fid = fopen(TEST.summaryFile,'a');                                     % Reopen the .csv summary
for i = 1:numel(target)
    [entryString,objectFlag] = scanObject(TEST,i,target{i});                % Get the report entry
    fprintf(TEST.fid,'%s\r\n',entryString);                                 % write the string to a file
    if objectFlag
        isIssue = 1;                                                        % Simply confirm if there are any issues
    end
end
fclose(TEST.fid);                                                           % Close the .csv summary
% /////////////////////////////////////////////////////////////////////////

% DELETE UNECESSARY DATA
[dirFlag,message] = rmdir(TEST.tempDir,'s');
if dirFlag
    fprintf('[%s]\tTemporary data removed @%s.\n',TEST.phase,TEST.tempDir);
else
    warning('Temporary still present at @%s.\n message:%s',TEST.outputPath,message);
end

% CONVERT THE .CSV TO A WORKBOOK
[~,~,a] = xlsread(TEST.summaryFile);
a = cellfun(@num2str, a,'UniformOutput',false);
xlswrite([strrep(TEST.summaryFile,'.csv',''),'.xlsx'],a);

% Delete the .csv
delete(TEST.summaryFile);
end

% Scan the complete directory (assumed filled with objects)
function [entryString,flag] = scanObject(TEST,i,objectLabel)

% DELCARE THE OPERATION DETAILS
fprintf('[%s]\tTest %d - object "%s"\n',TEST.phase,i,objectLabel);

% DEFAULT CONDITIONS
importStatus = 0;
simStatus = 0;
importMessage = '-';
simMessage = '-';
notes = '-';
flag = 0;

if contains(objectLabel,'.m')
    objectLabel = strrep(objectLabel,'.m',''); % Remove extension if necessary
end
    
% CREATE THE OBJECT SET
objectIndex = cell(TEST.objectNumber,1);
try
    for j = 1:TEST.objectNumber
        objectIndex{j} = eval(objectLabel);                                % Attempt to construct the object
    end
    importStatus = 1;
catch objectError
    warning(objectError.message);
    importMessage = objectError.message;
end

% ASSESS SIMULATION CONDITION
if ~importStatus
    fprintf('[%s]\tIMPORT FAILED: OBJECT(%s).\n',TEST.phase,objectLabel);
end

% IF THE OBJECT CAN BE SIMULATED, ATTEMPT SIMULATION
fprintf('[%s]\tTest %d - object "%s"\n',TEST.phase,i,objectLabel);
fprintf('[%s]\tSIMULATING...\n',TEST.phase);
try
    [~,~] = OMAS_initialise('objects',objectIndex,...
        'duration',1,...
        'dt',0.1,...
        'figures','None',...
        'verbosity',0,...
        'outputPath',TEST.tempDir);
    simStatus = 1;
catch simError
    warning(simError.message);
    simMessage = simError.message;
    notes = strcat('[function] -> ',simError.stack(1).name);
    notes = strcat(notes,' [line] -> ',num2str(simError.stack(1).line));
end

% Indicate if errors occurred
if importStatus || simStatus
   flag = 1; 
end

% CREATE SUMMARY ENTRY
if importStatus
    importStatus = 'PASSED';
else
    importStatus = 'FAILED';
end
if simStatus
    simStatus = 'PASSED';
else
    simStatus = 'FAILED';
end

% CREATE ROW ENTRY
rowData = cell(1,6);
rowData{1} = num2str(i);
rowData{2} = objectLabel;
rowData{3} = importStatus;
rowData{4} = importMessage;
rowData{5} = simStatus;
rowData{6} = simMessage;
rowData{7} = notes;
% Define entry to the log
entryString = rowData{1};
for j = 2:numel(rowData)
    stringValue = strrep(rowData{j},',','?');
    entryString = [entryString,',',stringValue];
end
end

% Get the XLS header
function [TEST] = getXlsHeader(TEST)
% PREPARE THE OVERVIEW FILE
rowData = {'#',...
           'name',...
           'import status',...
           'import message',...
           'OpenMAS status',...
           'OpenMAS message',...
           'notes'};
header_string = rowData{1};
for i = 2:length(rowData)
    header_string = [header_string,',',rowData{i}];
end
% Write the header to the file
fprintf(TEST.fid,'%s\r\n',header_string);
end