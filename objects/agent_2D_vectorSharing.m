%% GEOMETRIC COLLISION AVOIDANCE AGENT (agent_2D_vectorSharing.m) %%%%%%%%%
% This programs contains a generic agent object with the 2017 3D geometric
% avoidance alogorithm applied to conduct course corrections.

% Author: James A. Douthwaite

classdef agent_2D_vectorSharing < agent_2D & agent_vectorSharing
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods 
        % Constructor
        function [obj] = agent_2D_vectorSharing(varargin)

            % Call the super class
            obj@agent_2D(varargin);    

            % //////////////////// SENSOR PARAMETERS //////////////////////
%             [obj.SENSORS] = obj.GetDefaultSensorParameters();       % Default sensing
            [obj.SENSORS] = obj.GetCustomSensorParameters();       % Experimental sensing
            % /////////////////////////////////////////////////////////////

            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin); 
            % /////////////////////////////////////////////////////////////
        end
        % Setup - X = [x y psi dx dy dpsi]
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function calculates the intial state for a generic
            % object.
            % The default state vector:
            % [x y psi dx dy dpsi]
            % INITIALISE THE 2D STATE VECTOR WITH CONCANTINATED VELOCITIES
            [obj] = obj.initialise_2DVelocities(localXYZVelocity,localXYZrotations);
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            dt = ENV.dt;

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
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE
            % /////////////////////////////////////////////////////////////
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            % Design the current desired trajectory from the waypoint.
            desiredHeadingVector = obj.GetTargetHeading();
            desiredVelocity = desiredHeadingVector*obj.nominalSpeed;

            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  
            if ~isempty(avoidanceSet) 
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = obj.GetAvoidanceCorrection(desiredVelocity,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % ///////////////////// CONTROLLER ////////////////////////////
            [obj] = obj.controller(dt,desiredVelocity);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(ENV,algorithm_indicator,algorithm_dt);      % Record when the algorithm is ran
            obj.DATA.inputNames = {'$v_x$ (m/s)','$v_y$ (m/s)','$\dot{\psi}$ (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = obj.localState(4:6);         % Record the control inputs
            
            % // DISPLAY CONFLICT RESOLUTION
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                obj = obj.GetAnimationFrame(overHandle,ENV,'resolutionZone.gif');
                close(overHandle);
            end
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods (Access = public)
        % GET THE AVOIDANCE CORRECTION
        function [headingVector,speed] = GetAvoidanceCorrection(obj,desiredVelocity,visualiseProblem)
            % This function calculates the collision avoidance velocity in
            % light of the current obstacles
            
            % Check we aren't stopping
            if iszero(desiredVelocity)
               headingVector = [1;0];
               speed = 0;
               return
            end
            
            % AGENT KNOWLEDGE
            [p_i,v_i,r_i] = obj.GetAgentMeasurements(); % Its own position, velocity and radius
            
            % Define the obstacle list
            obstacleIDs = [obj.MEMORY([obj.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
                        
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            optimalSet = [];
            for j = 1:numel(obstacleIDs)
                % Fetch the obstacle trajectories
                p_j = obj.GetLastMeasurementByID(obstacleIDs(j),'position');
                
                % Neighbour conditions
                neighbourConditionA = j < obj.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < obj.neighbourDist;       % [CONFIRMED]
                if ~neighbourConditionA || ~neighbourConditionB
                    continue
                end
                
                % Get further obstacle information
                v_j = obj.GetLastMeasurementByID(obstacleIDs(j),'velocity');
                r_j = obj.GetLastMeasurementByID(obstacleIDs(j),'radius');
                
                % OBSTACLE KNOWLEDGE
                p_j = p_j + p_i;
                v_j = v_j + v_i;                                           % Convert relative parameters to absolute

                % COMPUTE THE VECTOR SHARING PROBLEM
                [Voptimal] = obj.define2DVectorSharingVelocity(desiredVelocity,p_i,v_i,r_i,p_j,v_j,r_j,visualiseProblem);
                optimalSet = horzcat(optimalSet,[Voptimal;abs(norm(desiredVelocity - Voptimal))]);
            end
            
            % INTERPRETING THE AVOIDANCE VELOCITY SET
            if ~isempty(optimalSet)
                inputDim = 2;
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
                headingVector = [1;0];   % Retain forward direction
                speed = 0;
            else
                speed = norm(avoidanceVelocity);
                headingVector = avoidanceVelocity/speed;
            end
        end
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
            % v_b - The obstacle's absolute velocity
            % r_b - The obstacle's radius
            % visualiseProblem - Plot flag
            % OUTPUTS:
            % U_a - The optimal vector heading, scaled by the desired
            %       velocity.
                        
            % Generate the 3D inputs 
            p_a = [p_a;0]; v_a = [v_a;0];
            p_b = [p_b;0]; v_b = [v_b;0];
            % Pass to the 3D function
            [U_a] = obj.define3DVectorSharingVelocity(desiredVelocity,...
                                                      p_a,v_a,r_a,p_b,v_b,r_b,...
                                                      visualiseProblem);
            % Reform the inputs for 2D application
            U_a = U_a(1:2,1);
        end
    end
    % The function uses the same sensor characteristics from the parent
    % class "agent_vectorSharing.m".
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]