%% TEST AGENT 
% This agent is used to examine the simulation feedback.

% Author: James A. Douthwaite

classdef agent_test_velocities < agent
    
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % CONTROL PARAMETERS
        priorError = zeros(3,1);
        % MISC TEST PARAMETERS
        controlVector;      % Properties for holding controller values
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_test_velocities(varargin)
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
            obj.sensorRange = 150;                                         % Virtual sensor range
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
            % UPDATE THE AGENT WITH THE NEW INFORMATION
            [obj,obstacleSet,~] = obj.getAgentUpdate(varargin{1});  

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            if ~isempty(obj.targetWaypoint)
                % DEFAULT THE HEADING VECTOR TOWARDS THE WAYPOINT
                waypointPosition = obj.targetWaypoint.state(1:3);
                desiredVelocity  = (waypointPosition/norm(waypointPosition))*norm(obj.localState(7:9));
            else
                desiredVelocity = obj.localState(7:9);
            end

            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
                if ~isempty(obstacleSet)
                    % PLOT THE OBSTACLE SCENE
                    [overHandle] = plotObstacleScene(obj,obstacleSet,overHandle);
                end
            end            

            %% \\\\\\\\\\\\\\ CALCULATE CONTROL INPUTS \\\\\\\\\\\\\\\\\\\\
            desiredHeading = desiredVelocity/norm(desiredVelocity);
            desiredSpeed = norm(desiredVelocity);
            [e_theta,e_lambda] = obj.getControlInputs([1;0;0],desiredHeading);
            e_control = [0;e_theta;-e_lambda]; % -ve in the NED control frame
            
            % BUILD THE CONTROL INPUTS
            Kp = 0.5; Kd = 0;  
            dV = Kp*eye(3)*e_control + Kd*eye(3)*(e_control-obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = e_control;
            
            %% \\\\\\\\\\\\\\ APPLY KINEMATIC LIMITS \\\\\\\\\\\\\\\\\\\\\\
            
            
            % \\\\\\\\ GET THE STATE UPDATE FROM THE CONTROL INPUTS \\\\\\\
            agentStateUpdate = obj.stateDynamics_velocities(dt,[desiredSpeed;0;0],dV);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                [agentStateUpdate] = obj.freezeAgent();
            end 
            % \\\\\\\\\\\\\\\ UPDATE THE CLASS GLOBAL PROPERTIES \\\\\\\\\\
            obj = obj.updateGlobalProperties(dt,agentStateUpdate);
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.inputs(1:3,TIME.currentStep) = dV;             % Record the control inputs 
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