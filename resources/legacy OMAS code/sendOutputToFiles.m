% SAVE DATA TO FILES
function sendOutputToFiles(META,EVENTS,DATA,objectIndex)
% This function is designed to handle the output data from the simulation
% and export the variables to background files
% INPUTS:
% META        - Local copy of the META variable
% EVENTS      - The cell array of event history objects
% DATA        - The output data structure
% objectIndex - The object (non-META) object class cell array

% DISPLAY ALL VARIABLES
whos;
% SAVE META DATA
save(strcat(META.outputPath,'META.mat'),'META');
save(strcat(META.outputPath,'OBJECTS.mat'),'objectIndex');
% SAVE EVENT HISTORY
save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
% SAVE OUTPUT DATA
save(strcat(META.outputPath,'DATA.mat'),'DATA');
% DISPLAY NOTIFICATION
fprintf('[%s]\tData objects outputted to file:\n',META.phase);
fprintf('[%s]\tDirectory: %s\n',META.phase,META.outputPath);
end

% % UPDATE OUTPUT DATA STRUCTURES
% META = META;                % Redefine for consistancy with the main cycle. 
% % SAVE META DATA
% save(strcat(META.outputPath,'META.mat'),'META');
% % SAVE THE OBJECT-INDEX
% save(strcat(META.outputPath,'OBJECTS.mat'),'objectIndex');
% % SAVE EVENT HISTORY
% save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
% % SAVE OUTPUT DATA
% save(strcat(META.outputPath,'DATA.mat'),'DATA');
% 
% clearvars META

% /////////////////// UPDATE THE SAVED DATA STRUCTURES ////////////////////
% UPDATE OUTPUT DATA STRUCTURES
% META = META;                % Redefine for consistancy with the main cycle. 
% % SAVE META DATA
% save(strcat(META.outputPath,'META.mat'),'META');
% % SAVE THE OBJECT-INDEX
% save(strcat(META.outputPath,'OBJECTS.mat'),'objectIndex');
% % SAVE EVENT HISTORY
% save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
% % SAVE OUTPUT DATA
% save(strcat(META.outputPath,'DATA.mat'),'DATA');
% clearvars META