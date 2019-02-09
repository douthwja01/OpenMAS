%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_2D_vectorSharing < agent_2D & agent_vectorSharing
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function [obj] = agent_2D_vectorSharing(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D(varargin);                                        
                                  
            % //////////////////// SENSOR PARAMETERS //////////////////////
            [obj] = obj.getDefaultSensorParameters();       % Default sensing
            %[obj] = obj.getCustomSensorParameters();       % Experimental sensing
            % /////////////////////////////////////////////////////////////

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % ///////////////////// SETUP FUNCTION ////////////////////////////
        % SETUP - X = [x y psi dx dy dpsi]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            % INITIALISE THE 2D STATE VECTOR WITH CONCANTINATED VELOCITIES
            [obj] = obj.initialise_2DVelocities(localXYZVelocity,localXYZrotations);
        end
        % ///////////////////// AGENT MAIN CYCLE //////////////////////////
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
            desiredHeadingVector = [1;0];
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
            
            % ///////////////////// CONTROLLER ////////////////////////////
            [obj] = obj.controller(dt,desiredVelocity);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputNames = {'Vx (m/s)','Vy (m/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = obj.localState(4:6);         % Record the control inputs
            
            % // DISPLAY CONFLICT RESOLUTION
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                obj = obj.getAnimationFrame(overHandle,TIME,'resolutionZone.gif');
                close(overHandle);
            end
        end
    end
    % ///////////////////// THE VECTOR SHARING METHODS ////////////////////
    methods (Access = public)
%          % GET THE AVOIDANCE CORRECTION
%         function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,knownObstacles,visualiseProblem)
%             % This function calculates the collision avoidance velocity in
%             % light of the current obstacles
%             
%             % AGENT KNOWLEDGE
%             [p_a,v_a,r_a] = obj.getAgentMeasurements(); % Its own position, velocity and radius
%            
%             % MOVE THROUGH THE PRIORITISED OBSTACLE SET 
%             optimalSet = []; 
%             for item = 1:length(knownObstacles)
%                 % NEIGHBOUR CONDITIONS
%                 neighbourConditionA = item < obj.maxNeighbours;                   % Maximum number of neighbours
%                 neighbourConditionB = norm(knownObstacles(item).position) < obj.neighbourDist;      % [CONFIRMED] 
%                 if ~neighbourConditionA || ~neighbourConditionB
%                     continue
%                 end
%                 % OBSTACLE KNOWLEDGE
%                 p_b = knownObstacles(item).position + p_a; 
%                 v_b = knownObstacles(item).velocity + v_a;                 % Convert relative parameters to absolute
%                 r_b = knownObstacles(item).radius;
%                 
%                 % COMPUTE THE VECTOR SHARING PROBLEM
%                 [Voptimal] = obj.defineOptimalVelocity(desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem);
%                 optimalSet = horzcat(optimalSet,[Voptimal;abs(norm(desiredVelocity - Voptimal))]);
%             end
%             
%             % INTERPRETING THE AVOIDANCE VELOCITY SET
%             inputDim = 2;
%             if ~isempty(optimalSet)
%                 % THE CLOSEST OBSTACLE
% %                 avoidanceVelocity = optimalSet(1:inputDim,1);
%                 % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
%                 [~,minIndex] = min(optimalSet((inputDim+1),:),[],2);       % Return the index of the smallest vector
%                 avoidanceVelocity = optimalSet(1:inputDim,minIndex);
%             else
%                 % NOTHING TO AVOID, CHOOSE OPTIMAL VELOCITY
%                 avoidanceVelocity = desiredVelocity;
%             end
% 
%             % CHECK VELOCITY IS PERMISSIBLE
%             if any(isnan(avoidanceVelocity))
%                 headingVector = [1;0];   % Retain forward direction
%                 speed = 0;
%             else
%                 speed = norm(avoidanceVelocity);
%                 headingVector = avoidanceVelocity/speed;
%             end
%         end
        % DEFINE THE VECTOR SHARING AVOIDANCE PROBLEM
        function [U_a] = define2DVectorSharingVelocity(obj,desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)
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
                r_m = r_m(1:2,1);
            end
            
            if norm(r_m) == 0 % No miss vector -> compensate
                r_m = 0.01*rand(2,1);
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
                q = quiver(p_a(1),p_a(2),v_a(1),v_a(2),'b'); q.AutoScaleFactor = 1;
                [Xa,Ya,Za] = obj.getSphere([p_a;0],r_a);
                sphZone = mesh(Xa,Ya,Za);
                set(sphZone,'facealpha',0.2,...
                            'FaceColor','b',...
                            'LineWidth',0.1,...
                            'EdgeAlpha',0.2);
                q = quiver(p_b(1),p_b(2),v_b(1),v_b(2),'r'); q.AutoScaleFactor = 1; 
                [Xb,Yb,Zb] = obj.getSphere([p_b;0],r_b);
                sphZone = mesh(Xb,Yb,Zb);
                set(sphZone,'facealpha',0.2,...
                            'FaceColor','r',...
                            'LineWidth',0.1,...
                            'EdgeAlpha',0.2);
                
                % FIGURE 3
                % PROJECTED POSITION AT tau ( NO INTERVENTION ) 
                q = quiver(p_a(1),p_a(2),VaDotTau(1),VaDotTau(2),'k'); q.AutoScaleFactor = 1;          % A's current trajectory at tau
                q = quiver(p_b(1),p_b(2),VbDotTau(1),VbDotTau(2),'k'); q.AutoScaleFactor = 1;          % A's current trajectory at tau
                
                % Near-miss vector @ tau
                q = quiver(vaTau(1),vaTau(2),r_m(1),r_m(2),'g'); q.AutoScaleFactor = 1;
                
                % SHARED VECTOR COMPONENTS
                q = quiver(vaTau(1),vaTau(2),r_vsa(1),r_vsa(2),'b'); q.AutoScaleFactor = 1;
                q = quiver(vbTau(1),vbTau(2),r_vsb(1),r_vsb(2),'r'); q.AutoScaleFactor = 1;
                               
                % THE OPTIMAL VELOCITIES
                q = quiver(p_a(1),p_a(2),U_a(1),U_a(2),'b'); q.AutoScaleFactor = 1;
                q = quiver(p_b(1),p_b(2),U_b(1),U_b(2),'r'); q.AutoScaleFactor = 1;
                               
                % DEFINE THE RESOLUTION ZONE
                % DRAW THE RESOLUTION ZONE
                scatter(conflictCenter(1),conflictCenter(2),'k');
                [X_zone,Y_zone,Z_zone] = obj.getSphere([conflictCenter;0],conflictRadius);
                sphZone = mesh(X_zone,Y_zone,Z_zone);
                set(sphZone,'facealpha',0.2,...
                            'FaceColor','g',...
                            'LineWidth',0.1,...
                            'EdgeAlpha',0.2);
                        
                % DRAW THE ACCELERATION
                q = quiver(Uatangent(1),Uatangent(2),r_m(1),r_m(2),'k'); q.AutoScaleFactor = 1;
                q = quiver(Ubtangent(1),Ubtangent(2),-r_m(1),-r_m(2),'k'); q.AutoScaleFactor = 1;
            end
        end
    end
    % The function uses the same sensor characteristics from the parent
    % class "agent_vectorSharing.m".
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]