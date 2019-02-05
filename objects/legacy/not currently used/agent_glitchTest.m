%% TEST AGENT 
% This agent is used to examine the simulation feedback.

% Author: James A. Douthwaite

classdef agent_glitchTest < agent
    properties
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_glitchTest(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'
            % PARAMETERISE KINEMATIC CONSTRAINTS
            obj.linearVelocityLimits = [30;30;1];                          % Limits on the agents velocity 
            obj.linearAccelerationLimits = [1;1;1];                        % Limits on the agents acceleration
            obj.angularVelocityLimits = [pi/3;pi/3;pi/3];
            obj.angularAccelerationLimits = [0.05;0.2;0.2]; % rad/s^2
            obj.sensorRange = inf;                                         % Virtual sensor range
            % DEFINE VIRTUAL PROPERTIES
            obj.VIRTUAL.detectionRange = obj.sensorRange;                  % Assign the range attribute to the SIM VIRTUAL attribute
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);

        end
        %% AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            observationSet = varargin{1}; % The detected objects
            
            % UPDATE THE AGENT WITH THE NEW INFORMATION
            [obj,~,~] = obj.getAgentUpdate(observationSet);  
           
            if ~isempty(obj.targetWaypoint)
            	obj.DATA.query(:,TIME.currentStep) = obj.targetWaypoint.state(1:3);
            else
                obj.DATA.query(:,TIME.currentStep) = zeros(3,1);
            end
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            if ~isempty(obj.targetWaypoint)
                % DEFAULT THE HEADING VECTOR TOWARDS THE WAYPOINT
                waypointPosition = obj.targetWaypoint.state(1:3);
                desiredVelocity  = (waypointPosition/norm(waypointPosition))*norm(obj.localState(7:9));
            else
                waypointPosition = [1;0;0];
                desiredVelocity = obj.localState(7:9);
            end
                        
            % CONTROLLER 
            [d_speed,d_heading] = obj.simpleController(desiredVelocity);       
            newHeadingRate = (d_heading)/dt;                               % The implied angular rates
            newSpeed = norm(obj.localState(7:9)) + d_speed;                % Absolute speed

            % SINGLE INTEGRATOR DYNAMICS
            newState = obj.stateDynamics_singleIntegrator(dt,newSpeed,newHeadingRate);
            
            % GLOBAL PROPERTIES UPDATE           
            obj = obj.updateGlobalProperties(dt,newState);        
            
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)','Roll Rate (rad/s)','Pitch Rate (rad/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [newSpeed;newState(4:6);newState(10:12)];         % Record the control inputs 
        end
    end
    %% PRIVATE METHODS
    methods
        % PLOT ALL THE OBSERVABLE OBSTACLES IN THE LOCAL FRAME
        function [figureHandle] = plotObstacleScene(obj,obstacleSet,figureHandle)
            
            % PLOT THE AGENT FOR PERSPECTIVE
            [X,Y,Z] = obj.getSphere(zeros(3,1),1);
            sphereHandle = mesh(X,Y,Z);
            set(sphereHandle,'facealpha',0.2,'FaceColor','g',...
                             'LineWidth',0.1,'EdgeAlpha',0.1); 
            % PLOT AGENT VELOCITY
            q = quiver3(0,0,0,obj.localState(4),obj.localState(5),obj.localState(6),'r');
            q.AutoScaleFactor = 1;
                            
            % MOVE THROUGH THE OBSTACLE SET AND PLOT THE OBSTACLE POSITIONS           
            for item = 1:length(obstacleSet)
               obstacle = obstacleSet(item);
               [X,Y,Z] = obj.getSphere(obstacleSet(item).state(1:3),obstacle.size);
               sphereHandle = mesh(X,Y,Z);
               set(sphereHandle,'facealpha',0.2,'FaceColor','r',...
                                'LineWidth',0.1,'EdgeAlpha',0.1); 
               % PLOT AGENT VELOCITY
               q = quiver3(obstacle.state(1),obstacle.state(2),obstacle.state(3)...
                          ,obstacle.state(4),obstacle.state(5),obstacle.state(6),'b');
               q.AutoScaleFactor = 1;
            end
        end
    end
    methods (Static)
        
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]