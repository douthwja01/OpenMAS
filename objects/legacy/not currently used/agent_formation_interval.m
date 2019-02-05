%% FORMATION - INTERVAL AVOIDANCE AGENT (agent_formation_interval.m) %%%%%%
% This agent utilises both the seperation base formation control agent and 
% the interval based collision avoidance agent to conduct coordinated
% manoeuvres.

% Author: James A. Douthwaite

classdef agent_formation_interval < agent_formation_seperation & agent_vectorSharing_interval
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR METHOD
        function obj = agent_formation_interval(varargin)           
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_formation_seperation(varargin);                      % Create the super class 'agent'            
            % CONTROL PARAMETERS
            obj.priorError = zeros(4,1);
            % PARAMETERISE KINEMATIC CONSTRAINTS
            obj.linearVelocityLimits = [2;2;2];                            % Limits on the agents velocity 
            obj.linearAccelerationLimits = [10;10;10];                     % Limits on the agents acceleration
            obj.angularVelocityLimits = [inf;inf;inf];
            obj.angularAccelerationLimits = [inf;inf;inf];
            % AVOIDANCE PARAMETERS  
            obj.obstacleSafetyFactor = 1.2;                                % Modify the apparent size of the obstacle               
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
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            observationSet = varargin{1};                                  % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet,~] = obj.getAgentUpdate(observationSet,dt);
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredSpeed = 0;
            desiredVelocity = [1;0;0]*desiredSpeed;
            
            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to define a desired trajectory using the formation 
            % controller that will cause the agents to converge on their
            % desired seperations.
            formationEnabled = 1; V = 0; 
            if ~isempty(agentSet) && formationEnabled
                % PASS AGENT SET TO FORMATION CONTROLLER
                [fi,V] = obj.formationControl_seperation(agentSet);        % Get the force vector
                % NORMALISE THE FORCE VECTOR
                norm_fi = norm(fi);
                unit_fi = fi/norm_fi;
                constrainedNorm = obj.boundValue(norm_fi,-obj.linearVelocityLimits(1),obj.linearVelocityLimits(1));
                % DEFINE VELOCITY REQUEST
                desiredVelocity = constrainedNorm*unit_fi;
            end                
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % We technically consider both obstacles and agents to be
            % obstacles. This collective set shall be referred to as the
            % 'avoidanceSet'
            avoidanceSet = vertcat(obstacleSet,agentSet);                  % Combine both the observed obstacles and agents in visual range
            % RUN THE AVOIDANCE ALGORITHM
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(desiredVelocity,avoidanceSet,visualiseProblem);
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % /////////////// PARSE CONTROLLER INPUTS /////////////////////
            desiredVelocity = mid(desiredVelocity);
            targetSpeed = norm(desiredVelocity);
            targetHeading = desiredVelocity/targetSpeed;                   % Get the control inputs
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            useController = 0;
            if useController
                [d_heading,d_speed,obj] = obj.trajectoryController(targetHeading,targetSpeed);  
            else
                [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);
            end
            newSpeed = (norm(obj.localState(7:9)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = d_heading;                                        % Heading change relative to straight ahead
            
            % //////////////// APPLY KINEMATIC LIMITS /////////////////////
            [newHeading,newSpeed] = obj.kinematicContraints(dt,newHeading,newSpeed);
            newVelocity = [1;0;0]*newSpeed;                                % Convert to local velocity vector
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            [newState] = obj.stateDynamics_velocities(dt,newVelocity,newHeading);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                newState = obj.freezeAgent();
            end
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:length(obj.priorError),TIME.currentStep) = [newSpeed;newState(4:6)];         % Record the control inputs 
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            obj.DATA.lypanov(TIME.currentStep) = V;                               % Record the lypanov value
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]