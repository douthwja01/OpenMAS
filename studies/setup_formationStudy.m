%% FORMATION CONTROL & COLLISION AVOIDANCE STUDY ( formationControl_VO.m) %
% This file conducts the study used to develope the formation
% control/collision avoidance algorithm used in the ___ journal article.

% Author: James A. Douthwaite

clear all; close all; 

% GET THE CURRENT DIRECTORY DETAILS
parentDir = 'studies';
cd ..
% ADD SIMULATION PATHS
addpath('environment');        
addpath('objects');  
addpath('scenarios'); 
addpath('Intlab_V7.1'); 

% LOAD THE INTERVAL ELEMENTS
wrkDir = pwd;
try 
    test = infsup(0,1);
    clearvars test 
    IntDir = strcat(pwd,'\Intlab_V7.1\startupJD.m');
catch 
    run('startup.m'); 
    cd(wrkDir)
end 

% //////////////////////////// BUILD AGENT SET ////////////////////////////
fprintf('|| Assigning agent definitions:\n');
agentNumber = 4;
for index = 1:agentNumber
    agentIndex{index} = agent_formation_VO();
end 

% ///////////////////// DEFINE THE FORMATION SCENARIO /////////////////////
ell_star = 5*double(~eye(agentNumber));
[ objectIndex ] = getScenario_formation('adjacencyMatrix',ell_star,...
                                        'agents',agentIndex,....
                                        'radius',10,...
                                        'velocities',0,...
                                        'repeatable',1,...                 % Random or repeatable
                                        'plot',1);

clearvars -except objectIndex

% /////////////////// DEFINE THE SIMULATION PARAMETERS ////////////////////
% figureSet = {'all'};
% % figureSet ={'eventoverview','seperations'};%'eventoverview',,'fig'};
figureSet = {'events','isometric'};

% //////////// INITIALISE THE SIMULATION WITH THE OBJECT INDEX ////////////
% [DATA,META] = simulation_initialise('objects',objectIndex,...
%                                     'simTime',15,...  
%                                          'dt',0.01,... 
%                                     'figures',figureSet,...
%                                  'threadPool',0);
%                              
% % ///////////////// STUDY SPECIFIC OUTPUT DATA HANDLING ///////////////////
% clearvars -except DATA META
% load(strcat(META.outputPath,'META.mat'));
% load(strcat(META.outputPath,'EVENTS.mat'));
% 
% % PLOT THE LYPANOV FUNCTION
% Vglobal = zeros(size(DATA.timeVector));
% for ind = 1:DATA.totalObjects
%     Vglobal = Vglobal + DATA.objectIndex{ind}.DATA.lypanov;
% end
% figureHandle = figure('Name','System Lypanov Function');
% plot(DATA.timeVector,Vglobal);
% annotationText = sprintf('Final value: %s',num2str(Vglobal(1,end)));
% text(DATA.timeVector(1,end),Vglobal(1,end),annotationText);
% savefig(figureHandle,strcat(META.outputPath,'Lypanov.fig')); 
% fprintf('..Complete.\n');