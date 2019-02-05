%% DEBUG SIMULATION SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation defines the simulation input conditions and agent
% scenario for the UKACC collision avoidance comparison examples.

clear all; close all;

% ADD THE PROGRAM PATHS
addpath('environment');
addpath('objects');  
addpath('scenarios'); 

% CODE DEBUG COMMANDS %
% profile on
% profile viewer

%If intlab needs to be reloaded
% addpath('toolboxes/Intlab_V7.1');   
% try 
%     wrkDir = pwd;
%     test = infsup(0,1);
%     clearvars test 
%     IntDir = strcat(pwd,'\Intlab_V7.1\startupJD.m');
% catch 
%     run('startup.m'); 
%     cd(wrkDir)
% end

fprintf('[SETUP]\tInitialising debug script.\n');

%% SIMULATION PARAMETERS
[~, userdir]   = system('echo %USERPROFILE%'); % Get desktop path
sim_outputPath = strcat(userdir,'\desktop\OpenMAS_data');
sim_vebosity   = 1;
sim_warningDistance = 2;
sim_maxDuration = 10; 
sim_timeStep   =  0.1;                 % RVO timestep (0.25s)

% sim_figureSet = {'all'};
sim_figureSet = {'fig','gif','inputs'};
% sim_figureSet = {'gif'}; 

%% SCENARIO PARAMETERS 
sim_agentNumber     = 5;                   % TOTAL NUMBER OF AGENTS
sim_agentOrbit      = 5; 
sim_waypointOrbit   = 10;
sim_offsetAngle     = 0;
sim_agentVelocity   = 0;
sim_obstacleNumber  = 4;
sim_noiseSigma      = 0.2;
sim_plotScenario    = logical(false);
sim_adjacencyMatrix = double(~eye(sim_agentNumber));

%% INITIALISE AGENTS
fprintf('[SETUP]\tAssigning agent definitions:\n');
for index = 1:sim_agentNumber
    % BASIC CLASSES
%     agentIndex{index} = objectDefinition();     
%     agentIndex{index} = agent();
%     agentIndex{index} = agent_2D_test();
%     agentIndex{index} = agent_example();

% DYNAMICS
%     agentIndex{index} = ARdrone_LQR();
%     agentIndex{index} = quadcopter();
%     agentIndex{index} = quadcopter_backup();
%     agentIndex{index} = quadcopter_shiyu();

%     agentIndex{index} = fixedWing();
%     agentIndex{index} = globalHawk();
%     agentIndex{index} = boeing737();

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
%     agentIndex{index} = agent_vectorSharing();
%     agentIndex{index} = agent_vectorSharing_interval();

% 2D COLLISION AVOIDANCE
%     agentIndex{index} = agent_2D();    
%     agentIndex{index} = agent_2D_VO();
%     agentIndex{index} = agent_2D_VO_withComplex();
%     agentIndex{index} = agent_2D_RVO();
%     agentIndex{index} = agent_2D_HRVO();
%     agentIndex{index} = agent_2D_RVO2(); 
%     agentIndex{index} = agent_2D_vectorSharing();
%     agentIndex{index} = agent_2D_vectorSharing_interval();

% AR-DRONE BASED COLLISION AVOIDANCE
%     agentIndex{index} = ARdrone_formation();

    % OTHER
%     agentIndex{index} = obstacle();
    agentIndex{index} = obstacle_cuboid();
%     agentIndex{index} = obstacle_spheroid();
%     agentIndex{index} = earth();

end

%% PLACE AGENT OBJECTS IN PRE-DEFINED SCENARIO
% FORMATION CONTROL TESTS
% [ objectIndex ] = getScenario_corridor('agents',agentIndex,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_formation_split('agents',agentIndex,'agentSpacing',4,'adjacencyMatrix',sim_adjacencyMatrix,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_formation_fourAgentsFourObstacles('agents',agentIndex,'obstacleRadius',2,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_fourCuboidObstacles('agents',agentIndex,'obstacleRadius',2,'plot',sim_plotScenario);

% OBSTACLE TESTS
% [ objectIndex ] = getScenario_fourCuboidObstacles('agents',agentIndex,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_obstacleTrack('agents',agentIndex,'plot',sim_plotScenario);

% AGENT TESTS 
[ objectIndex ] = getScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',sim_offsetAngle,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_concentricRing('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario,'noiseFactor',sim_noiseSigma);
% [ objectIndex ] = getScenario_concentricSphere('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_concentricAngle('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'angle',pi/4,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_SWIMSMART('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_UKACC('agents',agentIndex,'agentOrbit',sim_agentOrbit,'agentVelocity',sim_agentVelocity,'waypointOrbit',sim_waypointOrbit,'offsetAngle',pi/2,'plot',sim_plotScenario);

% WAYPOINT TESTS                            
% [ objectIndex ] = getScenario_waypointCurve('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);
% [ objectIndex ] = getScenario_waypoint_90Degrees('agents',agentIndex,'agentVelocity',sim_agentVelocity,'plot',sim_plotScenario);

% RL TESTS
% [ objectIndex ] = getScenario_earthOrbit();

%% %%%%%% INITIALISE THE SIMULATION WITH THE OBJECT INDEX %%%%%%%%%%%%%%%%%
[DATA,META] = OMAS_initialise('objects',objectIndex,...
                             'duration',sim_maxDuration,... 
                                   'dt',sim_timeStep,...
                              'figures',sim_figureSet,...
                      'warningDistance',sim_warningDistance,...
                            'verbosity',sim_vebosity,...
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
