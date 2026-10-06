%% 2D GEOMETRIC COLLISION AVOIDANCE AGENT (agent_2D_VO.m) %%%%%%%%%%%%%%%%%


% Author: James A. Douthwaite

classdef agent_2D_RVO < agent_2D_VO & agent_RVO
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE 2D RECIPROCAL VO METHOD
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTOR
        function obj = agent_2D_RVO(varargin)
            % Construct the agent object and initialise with the following
            % specific parameters.
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_VO(varargin);     
            
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [obj] = obj.ApplyUserOverrides(varargin); 
            % /////////////////////////////////////////////////////////////
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.
        
        % CALCULATE THE NECESSARY 2D AVOIDANCE VELOCITY
        function [headingVector,speed] = GetAvoidanceCorrection(obj,dt,desiredVelocity,visualiseProblem)
            % This function calculates the 2D avoidance velocity command
            % and returns it to be achieved by the controller.
            
            % Input sanity check
            assert(numel(desiredVelocity) == 2,'Desired velocity must be 2D');
            
            % AGENT KNOWLEDGE (2D)
            [p_i,v_i,r_i] = obj.GetAgentMeasurements();
            
            % Define the obstacle list
            obstacleIDs = [obj.MEMORY([obj.MEMORY.type] ~= OMAS_objectType.waypoint).objectID];
            
            
            % MOVE THROUGH THE PRIORITISED OBSTACLE SET
            VO = [];
            for item = 1:numel(obstacleIDs)
                % Get object data from memory structure
                p_j = obj.GetLastMeasurementByID(obstacleIDs(item),'position');
                v_j = obj.GetLastMeasurementByID(obstacleIDs(item),'velocity');
                r_j = obj.GetLastMeasurementByID(obstacleIDs(item),'radius'); 
                
                % Neighbour conditions
                neighbourConditionA = item < obj.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < obj.neighbourDist;       % [CONFIRMED] 
                neighbourConditionC = ~any(isnan(v_j));                    % Wait for a valid velocity reading
                if ~neighbourConditionB || ~neighbourConditionC
                    continue
                end

                % OBSTACLE KNOWLEDGE
                p_j = p_j + p_i; 
                v_j = v_j + v_i;                                           % Convert relative parameters to absolute
                tau_j = 0;
                
                % OBSTACLE TYPE BEHAVIOUR
                type_j = obj.GetLastMeasurementByID(obstacleIDs(item),'type');
                if OMAS_objectType.agent == type_j
                    % DEFINE RVO AS AGENT BEHAVIOUR
                    [VO_j] = obj.define2DReciprocalVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                else
                    % OBSTACLE BEHAVIOUR
                    [VO_j] = obj.define2DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                end       
                % CONCATINATE THE VO SET
                VO = [VO,VO_j];
            end            
            
            % THE CLEAR PATH STRATEGY
            [avoidanceVelocity] = obj.strategy_clearPath(v_i,desiredVelocity,VO,visualiseProblem);  

            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0];  % Retain previous heading
            end
            
            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem && obj.objectID == visualiseProblem
                % PLOT THE LEADING TANGENT VECTOR
                OMAS_axisTools.drawTriad([p_i;0],eye(3));
                % CURRENT VELOCITY
                q = quiver(gca,p_i(1),p_i(2),v_i(1),v_i(2),'m');
                q.AutoScaleFactor = 1;
                % DESIRED VELOCITY
                q = quiver(gca,p_i(1),p_i(2),desiredVelocity(1),desiredVelocity(2),'g');
                q.AutoScaleFactor = 1;   
                % AVOIDANCE VELOCITY
                q = quiver(gca,p_i(1),p_i(2),avoidanceVelocity(1),avoidanceVelocity(2),'b');
                q.AutoScaleFactor = 1;
            end
        end
    end
    % /////////////// RECIPRICOL VELOCITY OBSTACLE METHODS ////////////////
    methods (Access = public)
        % DEFINE THE 2D RECIPROCAL VELOCITY OBSTACLE (RVO)
        function [RVO] = define2DReciprocalVelocityObstacle(obj,p_i,v_i,r_i,p_j,v_j,r_j,tau_j,plotOn)
            % This function assembles the reciprocal velocity obstacle in 2D. 
            
            % MAP THE 2D INPUTS TO 3D 
            p_i = [p_i;0]; v_i = [v_i;0];
            p_j = [p_j;0]; v_j = [v_j;0];
            % CALL THE 3D VO GENERATION FUNCTION
            [RVO] = obj.define3DReciprocalVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,plotOn);
            % MAP THE 3D INPUTS BACK TO 2D INPUTS
            RVO.apex = RVO.apex(1:2,1);
            RVO.axisUnit = RVO.axisUnit(1:2,1);
            RVO.leadingEdgeUnit = RVO.leadingEdgeUnit(1:2,1);
            RVO.trailingEdgeUnit = RVO.trailingEdgeUnit(1:2,1);
        end
    end
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]