%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_vectorSharing.m) %%%%%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_vectorSharing < agent
%%% INITIALISE THE AGENT SPECIFIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % AVOIDANCE PARAMETERS
        neighbourDist = 15;  
        maxNeighbours = 10;
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function [this] = agent_vectorSharing(varargin)
            % Call the super class
            this@agent(varargin);                                    % Get super class 'agent'
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
            [this.SENSORS] = this.GetDefaultSensorParameters();       % Default sensing
            %[obj.SENSORS] = obj.GetCustomSensorParameters();       % Experimental sensing
            % /////////////////////////////////////////////////////////////
            
            % Assign defaults
            this.DYNAMICS = this.CreateDYNAMICS();                    % Dynamic
            this = this.SetBufferSize(1);                             % Only one needed
            this.detectionRadius = this.SENSORS.range;
            
            % //////////////// Check for user overrides ///////////////////            
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup - X = [x;x_dot]' 3D STATE VECTOR
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % BUILD THE STATE VECTOR FOR A 3D SYSTEM WITH CONCATINATED VELOCITIES
            [this] = this.setup_3DVelocities(localXYZVelocity,localXYZrotations);
        end        
        % Main
        function [this] = main(this,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
                        
            % PLOT AGENT FIGURE
            visualiseFlag = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseFlag == 1
                overHandle = figure(100+this.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
                        
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this] = this.GetAgentUpdate(ENV,varargin{1});                
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            headingVector   = this.GetTargetHeading();
            desiredVelocity = headingVector*this.v_nominal;
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            algorithm_start = tic; avoidanceEnabled = 1;
            if avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = this.GetAvoidanceCorrection(desiredVelocity,visualiseFlag);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
                       
            % /////// COMPUTE STATE CHANGE FROM CONTROL INPUTS ////////////
            [this] = this.Controller(ENV.dt,desiredVelocity);
                        
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            this.DATA.inputNames = {'$\phi$ (rad)','$\theta$ (rad)','$\psi$ (rad)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(4:6);  
        end
    end
    % AGENT SPECIFIC METHODS
    methods
        % GET THE AVOIDANCE CORRECTION
        function [heading,speed] = GetAvoidanceCorrection(this,desiredVelocity,visualiseProblem)
            % This function calculates the collision avoidance velocity in
            % light of the current obstacles
            
            % Check we aren't stopping
            [heading,speed] = this.nullVelocityCheck(desiredVelocity);
            if speed == 0
                return 
            end
            
            % AGENT KNOWLEDGE
            [p_i,v_i,r_i] = this.GetAgentMeasurements(); % Its own position, velocity and radius
            
            % Define the obstacle list
            obstacleIDs = [this.MEMORY([this.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
            
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            optimalSet = [];
            for item = 1:numel(obstacleIDs)
                % Get object data from memory structure
                p_j = this.GetLastMeasurementByID(obstacleIDs(item),'position');
                % Neighbour condition
                neighbourConditionA = item < this.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < this.neighbourDist;        
                if ~neighbourConditionA || ~neighbourConditionB
                    continue
                end
                
                % Get further obstacle information
                v_j = this.GetLastMeasurementByID(obstacleIDs(item),'velocity');
                r_j = this.GetLastMeasurementByID(obstacleIDs(item),'radius');
                
                % OBSTACLE KNOWLEDGE
                p_j = p_j + p_i;
                v_j = v_j + v_i; % Convert relative parameters to absolute

                % COMPUTE THE VECTOR SHARING PROBLEM
                [Voptimal] = this.Define3DVectorSharingVelocity(desiredVelocity,p_i,v_i,r_i,p_j,v_j,r_j,visualiseProblem);
                
                % Add the potential velocity and the control effort
                optimalSet = horzcat(optimalSet,[Voptimal;abs(norm(desiredVelocity - Voptimal))]);
            end
            
            % INTERPRETING THE AVOIDANCE VELOCITY SET
            if ~isempty(optimalSet)
                inputDim = 3;
                % FIND THE MINIMUM MAGNITUDE DEVIATION FROM THE DESIRED
                [~,minIndex] = min(optimalSet((inputDim+1),:),[],2);       % Return the index of the smallest vector
                avoidanceVelocity = optimalSet(1:inputDim,minIndex);
            else
                % NOTHING TO AVOID, CHOOSE OPTIMAL VELOCITY
                avoidanceVelocity = desiredVelocity;
            end
            
            % CHECK VELOCITY IS PERMISSIBLE
            if any(isnan(avoidanceVelocity))
                heading = [1;0;0];   % Retain forward direction
                speed = 0;
            else
                speed = norm(avoidanceVelocity);
                heading = avoidanceVelocity/speed;
            end
        end
        % DEFINE THE VECTOR SHARING AVOIDANCE PROBLEM
        function [U_a] = Define3DVectorSharingVelocity(this,desiredVelocity,p_a,v_a,r_a,p_b,v_b,r_b,visualiseProblem)
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
            
            % Define the current velocity as optimal
            U_a = desiredVelocity;
            
            % Input sanity #1 - No relative velocity.
            if ~any(v_b) || any(isnan(v_b))
                return
            end      
            
            % DEFINE THE INITIAL PROBLEM PARAMETERS
            r   = p_b - p_a;                          % The relative position vector = relative seperation
            c_b = v_b - v_a;                          % Define absolute velocity of b
            c_b_unit = c_b/norm(c_b);                 % Unit relative velocity of v_b
            
            % ///////////////// VECTOR SHARING PROBLEM ////////////////////
            
            % THE 'NEAR-MISS' VECTOR
            r_m = cross(c_b_unit,cross(r,c_b_unit));
            
            % CATCHA: No miss vector -> compensate
            if norm(r_m) == 0
                r_m = 0.01*rand(3,1);
            end
            
            % DEFINE THE TIME TO COLLISION (+ve converging)
            tau = -(dot(r,c_b)/dot(c_b,c_b));
            
            % COLLISION CHECK #2 - No future time to collision
            if tau < 0
                return
            end

            % VECTOR SHARING PROBLEM
            r_safe = r_a + r_b;         % Define the safe separation distance
            r_res = r_safe - norm(r_m); % Define the resolution zone
            
            % DEFINE THE VECTOR SHARING TERMS
            r_vsa = (norm(v_b)/(norm(v_a) + norm(v_b)))*(r_res/norm(r_m))*(-r_m);
            r_vsb = (norm(v_a)/(norm(v_a) + norm(v_b)))*(r_res/norm(r_m))*(r_m);
            
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_a = (v_a*tau + r_vsa);                                       % The ideal velocity to resolve collision
            U_b = (v_b*tau + r_vsb);
            
            % NORMALISED AND DEMENSION BY DESIRED SPEED
            U_a = (U_a/norm(U_a))*norm(desiredVelocity);
            
            % We can therefore write the optimal acceleration direction:
            % unitAcceleration = -r_m/norm(r_m);

            % PLOT THE SCENARIO IF REQUESTED
            if visualiseProblem == 1 && this.objectID == 1
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
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]