%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_2D_IA.m) %%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_2D_IA < agent_2D_vectorSharing & agent_interval 
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods
        % Constructor
        function [this] = agent_2D_IA(varargin)

            % CALL THE SUPERCLASS CONSTRUCTOR
            this@agent_2D_vectorSharing(varargin);                             % Get super class 'agent'
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
%             [obj.SENSORS] = obj.GetDefaultSensorParameters();       % Default sensing
            [this.SENSORS] = this.GetCustomSensorParameters();       % Experimental sensing
            % /////////////////////////////////////////////////////////////
            
            % Defaults
            this.v_nominal = 2;
            this.v_max = 4;
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
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
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+this.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            %fprintf('////////// LOOP START (%s) //////////\n',obj.name);
                       
            % /////////////// GET THE INFORMATION UPDATE //////////////////
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});        % IDEAL INFORMATION UPDATE
            % /////////////////////////////////////////////////////////////

            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            headingVector = this.GetTargetHeading();
            nominalVelocity = headingVector*this.v_nominal;
            % /////////////////////////////////////////////////////////////

            avoidanceVelocity = nominalVelocity;
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            algorithm_start = tic; algorithm_indicator = 0;
            if ~isempty([obstacleSet,agentSet])
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [avoidanceHeading,avoidanceSpeed] = this.GetAvoidanceCorrection(nominalVelocity,visualiseProblem);
                avoidanceVelocity = avoidanceHeading*avoidanceSpeed;
            end
            algorithm_dt = toc(algorithm_start);% Stop timing the algorithm
            % /////////////////////////////////////////////////////////////
                       
            if any(isnan(avoidanceVelocity))
                enactedVelocity = [0;0];
            else
                enactedVelocity = mid(avoidanceVelocity);
            end
                
            % ///////////// AGENT VELOCITY VECTOR CONTROLLER //////////////
            this = this.Controller(ENV.dt,enactedVelocity);

            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt); % Record when the algorithm is ran
            this.DATA.inputNames = {'$v_x (m/s)$','$v_y (m/s)$','$\dot{\psi} (rad/s)$'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(4:6);         % Record the control inputs
            
            %fprintf('////////// LOOP END (%s) //////////\n',obj.name);
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods    
       % GET THE INTERVAL SHARING VELOCITY CORRECTION
        function [headingVector,speed] = GetAvoidanceCorrection(obj,desiredVelocity,plotFlag)
            % This function is the highest level of the obstacle avoidance
            % algorithm. It accepts the agents current desired trajectory
            % and obstacles the agent is aware of (currently in memory) and
            % computes an alternative trajectory based on the given algorithm.
            
            % Sanity check #1 - Input
            assert(size(desiredVelocity,1) == 2,"Desired velocity must be a 2D vector.");
            if nargin < 3 
                plotFlag = 0;
            end
                        
            % Default behaviour 
            speed = norm(desiredVelocity);
            headingVector = desiredVelocity/speed;
                
            % Sanity check # 2 - Check we aren't stopping
            if sum(desiredVelocity) == 0
            	headingVector = [1;0]; % Ensure viable heading
            	return
            end
            
            % Parameters from agent i
            [p_i,v_i,r_i] = obj.GetAgentMeasurements();
            
            % Define the obstacle list
            obstacleIDs = [obj.MEMORY([obj.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
            
            U_candidates = intval([]);
            for i = 1:numel(obstacleIDs)
                % Fetch the obstacle trajectories
                p_j = obj.GetTrajectoryByID(obstacleIDs(i),'position');
                v_j = obj.GetTrajectoryByID(obstacleIDs(i),'velocity');
                r_j = obj.GetLastMeasurementByID(obstacleIDs(i),'radius');
                
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = i > obj.maxNeighbours;               % Maximum number of neighbours
                neighbourConditionB = inf(obj.inorm(p_j)) > obj.neighbourDist;   % [CONFIRMED] 
%                 neighbourConditionC = ~any(isnan(v_j));                 % Wait for a valid velocity reading
                if neighbourConditionA || neighbourConditionB 
                    continue
                end
                
%                 % For consistancy with other methods
%                 p_j = p_i + p_j;
%                 v_j = v_i + v_j;

                % Calculate the optimal avoidance routine
                U_i = obj.GetOptimalAvoidanceRegion(...
                    desiredVelocity,...
                    p_i,v_i,r_i,...
                    p_j,v_j,r_j,...
                    plotFlag);
                
                % Concatinate for storage
                U_candidates = horzcat(U_candidates,[U_i;obj.inorm(U_i)]);
            end
            
            % No executable velocities
            if size(U_candidates,2) < 1
                return
            end

            % Get the indices of the correction magnitudes in accending order
            [~,ind] = sort(mid(U_candidates(3,:)),2,'ascend'); 
            
            % Re-order the optimal regions prior to intersection
            U_candidates = U_candidates(1:2,ind);

            % Select the optimal region from the global set
            U_star = ones(2,1)*midrad(0,inf);
            for i = 1:size(U_candidates,2)
                % Intersect the current optimal with the next candidate
                temp = intersect(U_star,U_candidates(1:2,i));
                % The new common region is the optimal region
                if ~any(isnan(temp))
                    U_star = temp;
                end
            end
            
%             U_star(1) = mid(U_star(1));     % Assume the collision occurs early
%             U_star(2) = inf(U_star(2));         
            
            % Obtain control inputs from optimal region
            headingVector = obj.iunit(U_star);
            speed = obj.inorm(U_star);
        end
        % Get the optimal velocity region
        function [U_i] = GetOptimalAvoidanceRegion(obj,u_i,p_i,v_i,r_i,p_j,v_j,r_j,plotFlag)
            % This function computes the optimal avoidance region for
            % object j to be enacted by agent i.
            
            % INPUTS:
            % u_i - The desired velocity
            % p_i - Agent position
            % v_i - Agent velocity
            % r_i - Agent radius
            % p_j - Obstacle position
            % v_j - Obstacle velocity
            % r_j - Obstacle radii
            % OUTPUT:
            % U_i - The optimal avoidance region
            
            % Input sanity check
            assert(size(p_i,1) == size(p_j,1) & size(p_i,1) == 2,"Positions must be defined as [2x1] column vectors.");
            assert(size(v_i,1) == size(v_j,1) & size(v_j,1) == 2,"Velocities must be defined as [2x1] column vectors.");
            assert(size(r_i,1) == size(r_j,1) & size(r_j,1) == 1,"Radii must be defined as scalars.");
            
            % Define default behaviour
            U_i = u_i;          % Desired is initially the optimal
            % Relative properties
            r = p_j - p_i;      % Relative separation
            c = v_j - v_i;      % Relative velocity
            % Unit relative velocity
            c_unit = obj.iunit(c);
            
            % Define the time-to-collision (+ve converging)
            % tau = -(obj.idot(r,c)/obj.idot(c,c));
            A = obj.idot(r,c);
            B = obj.idot(c,c);
            tau = -(A/B);
%             tau = - obj.idivide(A,B); % tau = -(A/B);
            
            % The 'near-miss' vector
            r_m_unit = CWnormal_2D(c_unit);
            r_m_mag  = obj.idet(r,c_unit);
            r_m = r_m_mag*r_m_unit;
             
            % COLLISION CHECK #1 - Improper value
            if isnan(tau)
                warning("'tau' contains NaNs.");
                return
            end

            % COLLISION CHECK #2 - The maximal value indicates no chance of collision
            if sup(tau) < 0  
                fprintf("'tau' determined no chance of collision.\n");
                return
            end  
                        
            % Define the safe separation
            r_safe = r_i + r_j;                % Define the safe separation distance
            % Define the resolution region
            r_res  = r_safe - obj.inorm(r_m);  % Define the resolution zone
            
            % Calculate the share magnitude
            share = (obj.inorm(v_j)/(obj.inorm(v_i) + obj.inorm(v_j)));
            % Calculate the share projection
            r_vsi = share*r_res*(obj.iunit(-r_m));
            
            
            tau = inf(tau);           % Assume the collision occurs early
            r_vsi(2) = inf(r_vsi(2)); % Assume the maximal correction is needed

            
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_i = (v_i*tau + r_vsi);                                       % The ideal velocity to resolve collision
            
            % NORMALISED AND DEMENSION BY DESIRED SPEED
            speed   = obj.inorm(u_i);
            heading = obj.iunit(U_i);
            U_i = heading*speed;
        end
    end
end