%% COMPLETE OBJECT-SET CHECK FUNCTION (object_debug.m) %%%%%%%%%%%%%%%%%%%%
% This function attempts to insert an instance of all objects in 
clear all; close all;

% ADD THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('scenarios'); 

% TEST CONDITIONS
test_objectsPath = '\objects';
test_objectNumber = 2;
test_outputPath = [pwd,'\temp'];
test_phase = 'OBJECTS';
test_summaryFile = [pwd,'\objectSummary.csv'];

% PREPARE THE OVERVIEW FILE
rowData = {'Number','name','Import Status','Import Message','OpenMAS Status','OpenMAS Message'};
header_string = rowData{1};
for i = 2:length(rowData)
    header_string = [header_string,',',rowData{i}];
end

%write the string to a file
fid = fopen(test_summaryFile,'w');
fprintf(fid,'%s\r\n',header_string);
fclose(fid);

% GET ALL .m FILES IN THE OBJECTS DIRECTORIES
mFileList = dir([cd,test_objectsPath,'\*.m']);
% CREATE A TEST SCENARIO FOR EACH OBJECT IN THE SCENARIO
for i = 1:numel(mFileList)
    fileName = strrep(mFileList(i).name,'.m','');        % Get the object file name
    % DELCARE THE OPERATION DETAILS
    fprintf('[%s]\tTest %d - object "%s"\n',test_phase,i,fileName);
    % DEFAULT CONDITIONS    
    importStatus = 0;
    importMessage = '-';
    simStatus = 0; 
    simMessage = '-';
    
    % CREATE THE OBJECT SET
    objectIndex = cell(test_objectNumber,1);
    try
        for j = 1:test_objectNumber
            objectIndex{j} = eval(fileName);                 % Attempt to construct the object
        end
        importStatus = 1;
    catch objectError
        warning(objectError.message);
        importMessage = objectError.message;
    end
    
    % ASSESS SIMULATION CONDITION
    if ~importStatus
         fprintf('[%s]\tIMPORT FAILED: OBJECT(%s).\n',test_phase,fileName);
    else
        % IF THE OBJECT CAN BE SIMULATED, ATTEMPT SIMULATION
        fprintf('[%s]\tTest %d - object "%s"\n',test_phase,i,fileName);
        fprintf('[%s]\tSIMULATING...\n',test_phase);
        try
            [DATA,META] = OMAS_initialise('objects',objectIndex,...
                                         'duration',1,... 
                                               'dt',0.1,...
                                          'figures','None',...
                                        'verbosity',0,...
                                       'outputPath',test_outputPath);
            simStatus = 1; 
        catch simError
            warning(simError.message);
            simMessage = simError.message;
        end
    end
    
    %% CREATE SUMMARY ENTRY
    if importStatus == 1
        importStatus = 'PASSED';
    else
        importStatus = 'FAILED';
    end
    if simStatus == 1
        simStatus = 'PASSED';
    else
        simStatus = 'FAILED';
    end
    
    % CREATE ROW ENTRY
    rowData = cell(1,6);
    rowData{1} = num2str(i);
    rowData{2} = fileName;
    rowData{3} = importStatus;
    rowData{4} = importMessage;
    rowData{5} = simStatus;
    rowData{6} = simMessage;
    
    entryString = rowData{1};
    for j = 2:numel(rowData)
        stringValue = strrep(rowData{j},',','?');
        entryString = [entryString,',',stringValue];
    end
        
    %write the string to a file
    fid = fopen(test_summaryFile,'a');
    fprintf(fid,'%s\r\n',entryString);
    fclose(fid);
end

% DELETE UNECESSARY DATA
[flag,message] = rmdir(test_outputPath,'s');
if flag
    fprintf('[%s]\tTemporary data removed @%s.\n',test_phase,test_outputPath);
else
    warning('Temporary still present at @%s.\n message:%s',test_outputPath,message);
end

% CONVERT THE .CSV TO A WORKBOOK
[~,~,a] = xlsread(test_summaryFile);
a = cellfun(@num2str, a,'UniformOutput',false);
xlswrite([strrep(test_summaryFile,'.csv',''),'.xlsx'],a);
% REMOVE .CSV
delete(test_summaryFile); % Delete the original .csv