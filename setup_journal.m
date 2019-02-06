%% JOURNAL TOPIC SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program constructs the simulation example scenarios examined in the
% journal topic 

clear all; close all;

% ADD THE PROGRAM PATHS
addpath('scenarios'); 
addpath('environment');
addpath('objects');  

% CODE DEBUG COMMANDS %
% profile on
% profile viewer 

fprintf('[SETUP]\tInitialising the 2018 journal example scenarios.\n');

%% SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%');                             % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\OpenMAS_data');
sim_vebosity   = 1;
sim_warningDistance = 2;
sim_threadPool = logical(false);
sim_maxDuration = 30; 
sim_timeStep   =  0.25;                                                    % RVO timestep (0.25s)

% sim_figureSet = {'all'};
% sim_figureSet = {'plan','inputs','gif'};
sim_figureSet = {'gif','events'};

%% SCENARIO PARAMETERS 
sim_agentNumber = 1;                    % TOTAL NUMBER OF AGENTS
sim_agentOrbit  = 5; 
sim_waypointOrbit = 10;
sim_offsetAngle = 0;
sim_agentVelocity = 0;
sim_obstacleNumber = 4;
sim_noiseSigma = 0.2;
sim_plotScenario = 1;
sim_adjacencyMatrix = double(~eye(sim_agentNumber));

%% INITIALISE AGENTS
fprintf('[SETUP]\tAssigning agent definitions:\n');
for index = 1:sim_agentNumber
   
% 3D FORMATION CONTROL 
%     agentIndex{index} = agent_formation();
%     agentIndex{index} = agent_formation_boids();

%     agentIndex{index} = agent_formation_VO();
%     agentIndex{index} = agent_formation_RVO();
%     agentIndex{index} = agent_formation_HRVO();

% 2D FORMATION CONTROL
%     agentIndex{index} = agent_2D_formation_VO();
%     agentIndex{index} = agent_2D_formation_RVO();
%     agentIndex{index} = agent_2D_formation_HRVO();
%     agentIndex{index} = agent_2D_formation_RVO2();

% 3D COLLISION AVOIDANCE
%     agentIndex{index} = agent_VO();
%     agentIndex{index} = agent_RVO();
%     agentIndex{index} = agent_HRVO();

% 2D COLLISION AVOIDANCE
%     agentIndex{index} = agent_2D();    
%     agentIndex{index} = agent_2D_VO();
    agentIndex{index} = agent_2D_VO_withComplex();
%     agentIndex{index} = agent_2D_RVO();
%     agentIndex{index} = agent_2D_HRVO();
%     agentIndex{index} = agent_2D_RVO2(); 
end
% DESIGN AN OBSTACLE INDEX
% obstacleIndex = {};
% for index = 1:sim_obstacleNumber
%     obstacleIndex{index} = obstacle_cuboid();
% end

% FORMATION CONTROL TESTS
% [ objectIndex ] = getScenario_corridor('agents',agentIndex,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario);
[ objectIndex ] = getScenario_formation_split('agents',agentIndex,'agentSpacing',4,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_formation_fourAgentsFourObstacles('agents',agentIndex,'obstacleRadius',2,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_fourCuboidObstacles('agents',agentIndex,'obstacleRadius',2,'plot',sim_plotScenario);

% OBSTACLE TESTS
% [ objectIndex ] = getScenario_fourCuboidObstacles('agents',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_obstacleTrack('agents',agentIndex,'plot',sim_plotScenario);

% AGENT TESTS 
% [ objectIndex ] = getScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',sim_offsetAngle,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_concentricSphere('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_boidsExample('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypoints',10,'waypointOrbit',sim_waypointOrbit,'plot',sim_plotScenario);

%% %%%%%% INITIALISE THE SIMULATION WITH THE OBJECT INDEX %%%%%%%%%%%%%%%%%
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,... 
                                   'dt',sim_timeStep,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
                           'threadPool',sim_threadPool,...
                           'outputPath',sim_outputPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except DATA META
load(strcat(META.outputPath,'META.mat'));
load(strcat(META.outputPath,'EVENTS.mat'));

%% EXTRA GRAPHS (DEBUG)
% objectDATA = DATA.objectIndex{1,1}.DATA;
% figure()
% grid on; 
% plot(META.TIME.timeVector,objectDATA.query)
% xlabel('time (s)');
% legend('Location','southeast');
