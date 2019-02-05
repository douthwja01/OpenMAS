
% Initial Matlab script to test the initialisation procedure of the
% simulink simulation model and the modules imported from the system
% libraries.
% Author: James A. Douthwaite 20/10/2015

clear all; close all;
addpath('agents');    % Import vehicle class library

% Initialise N agents, cataloged in a managed manner, with a given state vector 
state = [1;2;3;4;5;6];
state2 = [7;8;9;10;11;12];
ID = zeros(10,2);
%ID = {'cat','dog'};
for i=1:length(ID)
    agents(i) = ARdrone('dog');%); %ID{i}             % Create a new agent class with ID
    
    agents(i) = agents(i).initialiseWithState(state);
    methods(agents(i))
    %     agents(i) = agents(i).initialise(state);
end

% display(agents); 
% pause;
% acceleration1 = [2;1;2];
% acceleration2 = [1;1;1];
% 
% for i = 1:length(ID)
%    display(agent(i).velocity);
%    agent(i) = agent(i).updateStateVector(acceleration1,acceleration2,1); 
%    display(agent(i).velocity); 
% end
