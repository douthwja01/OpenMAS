%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_IA.m) %%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_IA < agent_vectorSharing & agent_interval
    % This class contains the 'vectorSharing_interval' code that combines
    % aspects of the 'vectorSharing' algorithm with the principles of
    % interval analysis.
    
    properties
        U0 = ones(3,1)*midrad(0,inf);
    end
    
    %% ////////////////////////// MAIN METHODS ////////////////////////////
    methods (Access = public)
        % Constructor
        function [this] = agent_IA(varargin)
            % Call the super class
            this@agent_vectorSharing(varargin); 
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
            %[this.SENSORS] = this.GetDefaultSensorParameters();              % Default sensing
            [this.SENSORS] = this.GetCustomSensorParameters();                % Experimental sensing
            % /////////////////////////////////////////////////////////////
            
            % Kinematic limitations
            this.v_max = 4;
            this.v_nominal = 2;
            
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
            % this      - The updated project
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100+this.objectID);
                hold on; grid on;
                axis equal;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end
            
            %fprintf('////////// LOOP START (%s) //////////\n',this.name);
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});
            % /////////////////////////////////////////////////////////////
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            headingVector   = this.GetTargetHeading();
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
            
            enactedVelocity = mid(avoidanceVelocity);

            % ///////////// AGENT VELOCITY VECTOR CONTROLLER //////////////
            this = this.Controller(ENV.dt,enactedVelocity);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);         % Record the computation time
            this.DATA.inputNames = {...
                '$v_x$ (m/s)','$v_y$ (m/s)','$v_z$ (m/s)',...
                '$\dot{\phi}$ (rad/s)','$\dot{\theta}$ (rad/s)','$\dot{\psi}$ (rad/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(7:12); % Record the control inputs
        end
    end
    
    %% //////////////////////// AUXILLARY METHODS /////////////////////////
    methods
        % Avoidance routine (top level)
        function [headingVector,speed] = GetAvoidanceCorrection(this,desiredVelocity,plotFlag)
            % This function is the highest level of the obstacle avoidance
            % algorithm. It accepts the agents current desired trajectory
            % and obstacles the agent is aware of (currently in memory) and
            % computes an alternative trajectory based on the given algorithm.
            
            % Sanity check #1 - Input
            assert(size(desiredVelocity,1) == 3,"Desired velocity must be a 3D vector.");
            if nargin < 3
                plotFlag = 0;
            end
            
            % Default behaviour
            speed = norm(desiredVelocity);
            headingVector = desiredVelocity/speed;
            
            % Sanity check # 2 - Check we aren't stopping
            if sum(desiredVelocity) == 0
                headingVector = [1;0;0]; % Ensure viable heading
                return
            end
            
            % Parameters from agent i
            [p_i,v_i,r_i] = this.GetAgentMeasurements();
            
            % Define the obstacle list
            obstacleIDs = [this.MEMORY([this.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
            
            U_candidates = intval([]);
            for i = 1:numel(obstacleIDs)
                % Fetch the obstacle trajectories
                p_j = this.GetLastMeasurementByID(obstacleIDs(i),'position');

                % Neighbour conditions
                neighbourConditionA = i > this.maxNeighbours;                    % Maximum number of neighbours
                neighbourConditionB = inf(this.inorm(p_j)) > this.neighbourDist;  % Wait for a valid velocity reading
                if neighbourConditionA || neighbourConditionB
                    continue
                end
                
                % Get further obstacle information
                v_j = this.GetLastMeasurementByID(obstacleIDs(i),'velocity');
                r_j = this.GetLastMeasurementByID(obstacleIDs(i),'radius');
              
                % Calculate the optimal avoidance routine
                U_i = this.GetOptimalAvoidanceRegion(...
                    desiredVelocity,...
                    p_i,v_i,r_i,...
                    p_j,v_j,r_j);
                
                % Append the correction magnitude
                Unorm = this.inorm(U_i);
                % Concatinate for storage
                U_candidates = horzcat(U_candidates,[U_i;Unorm]);
            end
            
            % No executable velocities
            if size(U_candidates,2) < 1
                return
            end
            
            % Get the indices of the correction magnitudes in accending order
            [~,ind] = sort(mid(U_candidates(4,:)),2,'ascend');
            
            % Re-order the optimal regions prior to intersection
            U_candidates = U_candidates(1:3,ind);
            
            % Select the optimal region from the global set
            U_star = ones(3,1)*midrad(0,inf);
            for i = 1:size(U_candidates,2)
                % Intersect the current optimal with the next candidate
                temp = intersect(U_star,U_candidates(1:3,i));
                % The new common region is the optimal region
                if ~any(isnan(temp))
                    U_star = temp;
                end
            end
            
%             if mid(U_star(1)) > 0
%                 U_star(1) = inf(U_star(1));
%             else
%                 U_star(1) = sup(U_star(1));
%             end
            if mid(U_star(2)) > 0
                U_star(2) = inf(U_star(2));
            else
                U_star(2) = sup(U_star(2));
            end
            if mid(U_star(3)) > 0
                U_star(3) = inf(U_star(3));
            else
                U_star(3) = sup(U_star(3));
            end
            
            % Obtain control inputs from optimal region
            headingVector = this.iunit(U_star);
            speed = this.inorm(U_star);
        end
        % Get the optimal velocity region
        function [U_i] = GetOptimalAvoidanceRegion(this,u_i,p_i,v_i,r_i,p_j,v_j,r_j)
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
            assert(size(p_i,1) == size(p_j,1) & size(p_i,1) == 3,"Positions must be defined as [3x1] column vectors.");
            assert(size(v_i,1) == size(v_j,1) & size(v_j,1) == 3,"Velocities must be defined as [3x1] column vectors.");
            assert(size(r_i,1) == size(r_j,1) & size(r_j,1) == 1,"Radii must be defined as scalars.");
            
            % Define default behaviour
            U_i = u_i;          % Desired is initially the optimal
            % Relative properties
            r = p_j - p_i;      % Relative separation
            c = v_j - v_i;      % Relative velocity
            % Unit relative velocity
            c_unit = this.iunit(c);
            % The 'near-miss' vector
            r_m = this.icross(c_unit,this.icross(r,c_unit));
            
            % Define the time-to-collision (+ve converging)
            % tau = -(this.idot(r,c)/this.idot(c,c));
            A = this.idot(r,c);
            B = this.idot(c,c);
            tau = -(A/B);
%             tau = - this.idivide(A,B);
            
            % COLLISION CHECK #1 - Improper value
            if isnan(tau)
%                 warning("'tau' contains NaNs.");
                return
            end
            
            % COLLISION CHECK #2 - The maximal value indicates no chance of collision
            if sup(tau) < 0
%                 fprintf("'tau' determined no chance of collision.\n");
                return
            end
            
            % Define the safe separation
            r_safe = r_i + r_j;                     % Define the safe separation distance
            % Define the resolution region
            r_res  = sup(r_safe) - this.inorm(r_m);  % Define the resolution zone

            % Calculate the vector share region
%             r_vsa = (this.inorm(v_j)/(this.inorm(v_i) + this.inorm(v_j)))*(r_res/this.inorm(r_m))*(-r_m);

            % Calculate the share magnitude
            share = (this.inorm(v_j)/(this.inorm(v_i) + this.inorm(v_j)));
            % Calculate the share projection
            r_vsi = share*r_res*(-this.iunit(r_m));
            
            % Assume the worst
            tau = inf(tau);

            % Switching logic based on position of object
%             if mid(r_vsi(1)) > 0
%                 r_vsi(1) = inf(r_vsi(1));
%             else
%                 r_vsi(1) = sup(r_vsi(1));
%             end
%             if mid(r_vsi(2)) > 0
%                 r_vsi(2) = inf(r_vsi(2));
%             else
%                 r_vsi(2) = sup(r_vsi(2));
%             end
%             if mid(r_vsi(3)) > 0
%                 r_vsi(3) = inf(r_vsi(3));
%             else
%                 r_vsi(3) = sup(r_vsi(3));
%             end
            
            % CALCULATE THE CORRECTION UNIT VECTORS
            U_i = (v_i*tau + r_vsi);                                       % The ideal velocity to resolve collision
            
            % NORMALISED AND DEMENSION BY DESIRED SPEED
            speed   = this.inorm(u_i);
            heading = this.iunit(U_i);
            U_i = heading*speed;
        end
    end
end