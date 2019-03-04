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
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_VO(varargin);     
                        
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.
        
        % CALCULATE THE NECESSARY 2D AVOIDANCE VELOCITY
        function [headingVector,speed] = GetAvoidanceCorrection(obj,dt,desiredVelocity,knownObstacles,visualiseProblem)
            % This function calculates the 2D avoidance velocity command
            % and returns it to be achieved by the controller.
            
            % AGENT KNOWLEDGE
            p_a = obj.localState(1:2,1);
            v_a = obj.localState(4:5,1);
            r_a = obj.radius;
%             [p_a,v_a,r_a] = obj.GetAgentMeasurements();
            


            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:length(knownObstacles)
                % GENERAL NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;                                 % Maximum number of neighbours
                neighbourConditionB = norm(knownObstacles(item).position) < obj.neighbourDist;  % [CONFIRMED] 
                neighbourConditionC = ~any(isnan(knownObstacles(item).velocity));               % Wait for a valid velocity reading
                if ~neighbourConditionB || ~neighbourConditionC
                    continue
                end
                
                % GENERAL OBSTACLE PARAMETERS
                p_b = knownObstacles(item).position + p_a; 
                v_b = knownObstacles(item).velocity + v_a;                 % Convert relative parameters to absolute
                r_b = knownObstacles(item).radius;
                tau_b = 0;

                % OBSTACLE TYPE BEHAVIOUR
                if knownObstacles(item).type == OMAS_objectType.agent
                    % DEFINE RVO AS AGENT BEHAVIOUR
                    [VO_i] = obj.define2DReciprocalVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,visualiseProblem);
                else
                    % OBSTACLE BEHAVIOUR
                    [VO_i] = obj.define2DVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,visualiseProblem);
                end       
                % CONCATINATE THE VO SET
                VO = [VO,VO_i];
            end
            
            % GET THE CAPABLE VELOCITIES
%             capableVelocities = obj.feasabliltyMatrix;                     % Get all the feasible velocities this timestep 
            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
%             escapeVelocities = obj.GetEscapeVelocities(capableVelocities,VO);      % The viable velocities from the capable        
            % SEARCH THE VIABLE ESCAPE VELOCITIES FOR THE VELOCITY WITH THE 
            % SMALLEST DEVIATION FROM THE DESIRED VELOCITY
%             [avoidanceVelocity] = obj.strategy_minimumDifference(desiredVelocity,escapeVelocities);
            
            % THE CLEAR PATH STRATEGY
            [avoidanceVelocity] = obj.strategy_clearPath(v_a,desiredVelocity,VO,visualiseProblem);  

            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0];  % Retain previous heading
            end
            
            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem && obj.objectID == visualiseProblem
                % PLOT THE LEADING TANGENT VECTOR
                OMAS_axisTools.drawTriad([p_a;0],eye(3));
                % CURRENT VELOCITY
                q = quiver(gca,p_a(1),p_a(2),v_a(1),v_a(2),'m');
                q.AutoScaleFactor = 1;
                % DESIRED VELOCITY
                q = quiver(gca,p_a(1),p_a(2),desiredVelocity(1),desiredVelocity(2),'g');
                q.AutoScaleFactor = 1;   
                % AVOIDANCE VELOCITY
                q = quiver(gca,p_a(1),p_a(2),avoidanceVelocity(1),avoidanceVelocity(2),'b');
                q.AutoScaleFactor = 1;
            end
        end
    end
    % /////////////// RECIPRICOL VELOCITY OBSTACLE METHODS ////////////////
    methods (Access = public)
        % DEFINE THE 2D RECIPROCAL VELOCITY OBSTACLE (RVO)
        function [RVO] = define2DReciprocalVelocityObstacle(obj,p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn)
            % This function assembles the reciprocal velocity obstacle in 2D. 
            
            % MAP THE 2D INPUTS TO 3D 
            p_a = [p_a;0]; v_a = [v_a;0];
            p_b = [p_b;0]; v_b = [v_b;0];
            % CALL THE 3D VO GENERATION FUNCTION
            [RVO] = obj.define3DReciprocalVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn);
            % MAP THE 3D INPUTS BACK TO 2D INPUTS
            RVO.apex = RVO.apex(1:2,1);
            RVO.axisUnit = RVO.axisUnit(1:2,1);
            RVO.leadingEdgeUnit = RVO.leadingEdgeUnit(1:2,1);
            RVO.trailingEdgeUnit = RVO.trailingEdgeUnit(1:2,1);
        end
    end
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]