%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_vectorSharing < agent
%%% INITIALISE THE AGENT SPECIFIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % AGENT PARMETERS
        radius = 0.5;
        % AVOIDANCE PARAMETERS
        neighbourDist = 15;  
        maxNeighbours = 10;
    end
%%  CLASS METHODS
    methods
        % CONSTRUCTOR
        function obj = agent_vectorSharing(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'
            
            % DYNAMIC PARAMETERS
            [obj] = obj.getDynamicParameters();
            
            % SENSOR PARAMETERS
            [obj] = obj.getDefaultSensorParameters();
            %[obj] = obj.getCustomSensorParameters();
            
            % RE-ESTABLISH THE (SIMULATOR PARAMETERS)
            obj.VIRTUAL.radius = obj.radius; 
            obj.VIRTUAL.detectionRadius = obj.SENSORS.range;
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
        function [obj] = main(obj,TIME,varargin)
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
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            % /////////////////// DEFAULT BEHAVIOUR ///////////////////////
            desiredSpeed = obj.nominalSpeed;
            desiredHeadingVector = [1;0;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});       % IDEAL INFORMATION UPDATE
%             [obj,obstacleSet,agentSet] = obj.getSensorUpdate(dt,varargin{1}); % REALISTIC INFORMATION UPDATE
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            if ~isempty(obj.targetWaypoint)
                desiredHeadingVector = obj.targetWaypoint.position/norm(obj.targetWaypoint.position);
                desiredVelocity = desiredHeadingVector*desiredSpeed; % Desired relative velocity
            end
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                 [desiredHeadingVector,desiredSpeed] = obj.getAvoidanceCorrection(desiredVelocity,avoidanceSet,visualiseProblem);
                 desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
                       
            % /////// COMPUTE STATE CHANGE FROM CONTROL INPUTS ////////////
            [obj] = obj.controller(dt,desiredVelocity);
                        
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(TIME,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'Vx (m/s)','Pitch (rad)','Roll (rad)','Yaw (rad)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [obj.localState(7);obj.localState(4:6)];         % Record the control inputs 
        end
        % INITIALISE SENSORS (NOISY SENSOR PARAMETERS)
        function [obj] = getCustomSensorParameters(obj)
            % This function is designed to populate the SENSOR field with
            % representative sensor uncertainty.
            % ADD ADDITIONAL PARAMETERS TO THE SENSOR DESCRIPTIONS
            obj.SENSORS.range = inf;                                       % Assume the agent has perfect environmental knowledge (m)
            obj.SENSORS.sigma_position = 0.5;                              % Accurate to within 0.5m
            obj.SENSORS.sigma_velocity = 0.1;                              % Accurate to within 0.1m/s
            obj.SENSORS.sigma_rangeFinder = 0.1;                           % Accurate to within 0.1m
            obj.SENSORS.sigma_camera = 5.21E-5;                            % One pixel in a 1080p image
            obj.SENSORS.sampleFrequency = inf;                             % Sensing has perfect precision
        end
    end
    % AGENT SPECIFIC METHODS
    methods
        % GET THE AVOIDANCE CORRECTION
        function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,knownObstacles,visualiseProblem)
            % This function calculates the collision avoidance velocity in
            % light of the current obstacles
            
            % AGENT KNOWLEDGE
            [p_a,v_a,r_a] = obj.getAgentMeasurements(); % Its own position, velocity and radius
            
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            optimalSet = [];
            for item = 1:length(knownObstacles)
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;                   % Maximum number of neighbours
                neighbourConditionB = norm(knownObstacles(item).position) < obj.neighbourDist;      % [CONFIRMED] 
                if ~neighbourConditionA || ~neighbourConditionB
                    continue
                end
                % OBSTACLE KNOWLEDGE
                p_b = knownObstacles(item).position + p_a;
                v_b = knownObstacles(item).velocity + v_a;                 % Convert relative parameters to absolute
                r_b = knownObstacles(item).radius;
                
                % COMPUTE THE VECTOR SHARING PROBLEM
                [Voptimal] = obj.defineOptimalVelocity(desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem);
                optimalSet = horzcat(optimalSet,[Voptimal;abs(norm(desiredVelocity - Voptimal))]);
            end
            
            % INTERPRETING THE AVOIDANCE VELOCITY SET
            inputDim = 3;
            if ~isempty(optimalSet)
                % THE CLOSEST OBSTACLE
                %                 avoidanceVelocity = optimalSet(1:inputDim,1);
                % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
                [~,minIndex] = min(optimalSet((inputDim+1),:),[],2);       % Return the index of the smallest vector
                avoidanceVelocity = optimalSet(1:inputDim,minIndex);
            else
                % NOTHING TO AVOID, CHOOSE OPTIMAL VELOCITY
                avoidanceVelocity = desiredVelocity;
            end
            
            % CHECK VELOCITY IS PERMISSIBLE
            if any(isnan(avoidanceVelocity))
                headingVector = [1;0;0];   % Retain forward direction
                speed = 0;
            else
                speed = norm(avoidanceVelocity);
                headingVector = avoidanceVelocity/speed;
            end
        end
        % DEFINE THE VECTOR SHARING AVOIDANCE PROBLEM
        function [U_a] = defineOptimalVelocity(obj,desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)
            % This function calculates the avoidance vectors based on the
            % principle of vector sharing.
            % INPUTS:
            % desiredVelocity - The true optimal vector
            % p_a - The agent's absolute position
            % v_a - The agent's absolute velocity
            % r_a - The agent's radius
            % p_b - The obstacle's absolute position
            % p_b - The obstacle's absolute velocity
            % p_b - The obstacle's radius
            % visualiseProblem - Plot flag
            % OUTPUTS:
            % U_a - The optimal vector heading, scaled by the desired
            %       velocity.
            
            % DEFINE THE INITIAL PROBLEM PARAMETERS
            r   = p_b - p_a;                          % The relative position vector = relative seperation
            c_b = v_b - v_a;                          % Define absolute velocity of b
            c_b_unit = c_b/norm(c_b);                 % Unit relative velocity of v_b
                        
            % THE 'NEAR-MISS' VECTOR
            if numel(desiredVelocity) == 3
                r_m = cross(c_b_unit,cross(r,c_b_unit));
            else
                c_temp = [c_b_unit;0]; r_temp = [r;0];
                r_m = cross(c_temp,cross(r_temp,c_temp));
            end
            
            if norm(r_m) == 0 % No miss vector -> compensate
                r_m = 0.01*rand(3,1);
            end    
            
            % DEFINE THE TIME TO COLLISION (+ve converging)
            tau = -(dot(r,c_b)/dot(c_b,c_b));
            
            % AVOIDANCE CONDITION
            if tau < 0
                U_a = desiredVelocity;
                return
            end
            
            % DEFINE THE REST REGION
            r_safe = r_a + r_b; %;+ 0.5;
            
            % RESOLUTION ZONE
            r_res = r_safe - norm(r_m);                                    % Define the rest region
            
            % DEFINE THE VECTOR SHARING TERMS
            r_vsa = (norm(v_b)/(norm(v_a) + norm(v_b)))*(r_res/norm(r_m))*(-r_m);
            r_vsb = (norm(v_a)/(norm(v_a) + norm(v_b)))*(r_res/norm(r_m))*(r_m);
            
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (v_a*tau + r_vsa);                                       % The ideal velocity to resolve collision
            U_b = (v_b*tau + r_vsb);
            
            % NORMALISED AND DEMENSION BY DESIRED SPEED
            U_a = (U_a/norm(U_a))*norm(desiredVelocity);
            
            % DEFINE THE OPTIMAL ACCELERATION VECTOR
            % Based on the hamiltonian, H, the vector that minimises the
            % cost expression J, is an ahat that is parallel to -r_m:
            % J = -0.5*norm(r_m)^2;
            % H = dot(-r_m,c_b) + (a_mod*tau)*dot(-r_m,unit_acc);
            
            % We can therefore write the optimal acceleration direction:
            %             unitAcceleration = -r_m/norm(r_m);
            
            
            % PLOT THE SCENARIO IF REQUESTED
            if visualiseProblem == 1 && obj.objectID == 1
                % A & B's current trajectory at tau
                VaDotTau = v_a*tau;
                VbDotTau = v_b*tau;     % Velocity projections
                vaTau = VaDotTau + p_a;
                vbTau = VbDotTau + p_b; % The vector sharing parameters
                conflictCenter = vaTau + 1/2*r_m;
                conflictRadius = 0.5*(norm(r_m)) + 0.5*(norm(r_vsa) + norm(r_vsb));
                Uatangent = U_a + p_a;
                Ubtangent = U_b + p_b;
                
                hold on; grid on;
                axis equal;
                
                % PLOT THE TWO AGENTS & VELOCITIES
                q = quiver3(p_a(1),p_a(2),p_a(3),v_a(1),v_a(2),v_a(3),'b'); q.AutoScaleFactor = 1;
                [Xa,Ya,Za] = OMAS_axisTools.defineSphere(p_a,r_a);
                sphZone = mesh(Xa,Ya,Za);
                set(sphZone,'facealpha',0.2,...
                    'FaceColor','b',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0.2);
                q = quiver3(p_b(1),p_b(2),p_b(3),v_b(1),v_b(2),v_b(3),'r'); q.AutoScaleFactor = 1;
                [Xb,Yb,Zb] = OMAS_axisTools.defineSphere(p_b,r_b);
                sphZone = mesh(Xb,Yb,Zb);
                set(sphZone,'facealpha',0.2,...
                    'FaceColor','r',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0.2);
                
                % FIGURE 3
                % PROJECTED POSITION AT tau ( NO INTERVENTION )
                q = quiver3(p_a(1),p_a(2),p_a(3),VaDotTau(1),VaDotTau(2),VaDotTau(3),'k'); q.AutoScaleFactor = 1;          % A's current trajectory at tau
                q = quiver3(p_b(1),p_b(2),p_b(3),VbDotTau(1),VbDotTau(2),VbDotTau(3),'k'); q.AutoScaleFactor = 1;          % A's current trajectory at tau
                
                % Near-miss vector @ tau
                q = quiver3(vaTau(1),vaTau(2),vaTau(3),r_m(1),r_m(2),r_m(3),'g'); q.AutoScaleFactor = 1;
                
                % SHARED VECTOR COMPONENTS
                q = quiver3(vaTau(1),vaTau(2),vaTau(3),r_vsa(1),r_vsa(2),r_vsa(3),'b'); q.AutoScaleFactor = 1;
                q = quiver3(vbTau(1),vbTau(2),vbTau(3),r_vsb(1),r_vsb(2),r_vsb(3),'r'); q.AutoScaleFactor = 1;
                
                % THE OPTIMAL VELOCITIES
                q = quiver3(p_a(1),p_a(2),p_a(3),U_a(1),U_a(2),U_a(3),'b'); q.AutoScaleFactor = 1;
                q = quiver3(p_b(1),p_b(2),p_b(3),U_b(1),U_b(2),U_b(3),'r'); q.AutoScaleFactor = 1;
                
                % DEFINE THE RESOLUTION ZONE
                % DRAW THE RESOLUTION ZONE
                scatter3(conflictCenter(1),conflictCenter(2),conflictCenter(3),'k');
                [X_zone,Y_zone,Z_zone] = OMAS_axisTools.defineSphere(conflictCenter,conflictRadius);
                sphZone = mesh(X_zone,Y_zone,Z_zone);
                set(sphZone,'facealpha',0.2,...
                    'FaceColor','g',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0.2);
                
                % DRAW THE ACCELERATION
                q = quiver3(Uatangent(1),Uatangent(2),Uatangent(3),r_m(1),r_m(2),r_m(3),'k'); q.AutoScaleFactor = 1;
                q = quiver3(Ubtangent(1),Ubtangent(2),Ubtangent(3),-r_m(1),-r_m(2),-r_m(3),'k'); q.AutoScaleFactor = 1;
            end
        end
    end
    % //////////////////// CONTROLLER/ PLANT OPERATIONS ///////////////////
    methods
        % INITIALISE THE [x;x_dot]' 3D STATE VECTOR
        function [obj] = initialise_localState(obj,localXYZVelocity,localXYZrotations)
            % The state initialiser must be called 'initialise_localState'
            % and instead calls the 'initialise_3DVelocities' function in
            % this case. 
            
            % BUILD THE STATE VECTOR FOR A 3D SYSTEM WITH CONCATINATED VELOCITIES
            [obj] = obj.initialise_3DVelocities(localXYZVelocity,localXYZrotations);
        end
    end
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]