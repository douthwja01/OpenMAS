%% 2D GEOMETRIC COLLISION AVOIDANCE AGENT (agent_2D_HRVO.m) %%%%%%%%%%%%%%%


% Author: James A. Douthwaite

classdef agent_2D_HRVO < agent_2D_RVO & agent_HRVO 
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE HRVO METHOD
    end
%%  CLASS METHODS
    methods
        % CONSTRUCTOR
        function obj = agent_2D_HRVO(varargin)
            % Construct the agent object and initialise with the following
            % specific parameters.
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_RVO(varargin); 
                        
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.
        
        % CALCULATE THE NECESSARY 2D AVOIDANCE VELOCITY
        function [headingVector,speed] = getAvoidanceCorrection(obj,dt,desiredVelocity,knownObstacles,visualiseProblem)
            % This function calculates the 2D avoidance velocity command
            % and returns it to be achieved by the controller.
            
            % AGENT KNOWLEDGE
            p_a = obj.localState(1:2,1);
            v_a = obj.localState(4:5,1);
            r_a = obj.VIRTUAL.radius;
%             [p_a,v_a,r_a] = obj.getAgentMeasurements();
            
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
                    [VO_i] = obj.define2DHybridReciprocalVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,visualiseProblem);
                else
                    % OBSTACLE BEHAVIOUR
                    [VO_i] = obj.define2DVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,visualiseProblem);
                end
                % CONCATINATE THE VO SET
                VO = [VO,VO_i];
            end
            
%             % GET THE CAPABLE VELOCITIES
%             capableVelocities = obj.feasabliltyMatrix;
%             % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
%             escapeVelocities = obj.getEscapeVelocities(capableVelocities,VO);              % The viable velocities from the capable 
%             % APPLY THE MINIMUM DIFFERENCE SEARCH STRATEGY
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
                
%                 plotCapableVelocities = capableVelocities + p_a;
%                 plotEscapeVelocities = escapeVelocities + p_a;
%                 plotAvoidanceVelocities = avoidanceVelocity + p_a;
%                 scatter(gca,plotCapableVelocities(1,:),plotCapableVelocities(2,:),'r'); 
%                 scatter(gca,plotEscapeVelocities(1,:),plotEscapeVelocities(2,:),'g');                 
%                 scatter(gca,plotAvoidanceVelocities(1),plotAvoidanceVelocities(2),'b','filled'); 
            end
        end
    end    
    % /////////// HYBRID RECIPRICOL VELOCITY OBSTACLE METHODS /////////////
    methods (Access = public)
        % DEFINE THE 2D HYBRID RECIPROCAL VELOCITY OBSTACLE (HRVO)
        function [HRVO] = define2DHybridReciprocalVelocityObstacle(obj,p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn)
            % This function assembles the reciprocal velocity obstacle in 2D.

            % GENERATE A STANDARD VELOCITY OBSTACLE
            [VO]  = obj.define2DVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn);
            [RVO] = obj.define2DReciprocalVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn);
            
            % CALCULATE THE VECTOR INTERSECTIONS DEPENDING ON ORIENTATION
            % OF Va WITH RESPECT TO THE RVO CENTER-LINE
            if RVO.isVaLeading
                % Va is on the leading side of the RVO
                %apx is the intersection of the learding side of the RVO
                %and the trailing side of the VO.
                % INTERSECT LEADING RVO AND TRAILING VO
                RVOelement = RVO.leadingEdgeUnit;
                VOelement  = VO.trailingEdgeUnit;
                % DELARE NEW HRVO PROPERTES
                unit_leadingTangent  = RVOelement;
                unit_trailingTangent = VOelement;
            else
                % INTERSECT TRAILING RVO AND LEADING VO
                VOelement  = VO.leadingEdgeUnit;
                RVOelement = RVO.trailingEdgeUnit;
                % DELARE NEW HRVO PROPERTES
                unit_leadingTangent  = VOelement;
                unit_trailingTangent = RVOelement;
            end
            % GET THE INTERSECTION POINT
            [HRVOapex, isSuccessful] = obj.findAny2DIntersection(VO.apex,VOelement,RVO.apex,-RVOelement); % Project the RVO vector towards the VO element
            if ~isSuccessful
                HRVOapex = RVO.apex;
                warning('no line intersection');
            end
            % //////// CONSTRUCT THE HRVO FROM THE RVO TEMPLATE ///////////
            HRVO = RVO;
            HRVO.apex = HRVOapex;
            HRVO.leadingEdgeUnit = unit_leadingTangent;
            HRVO.trailingEdgeUnit = unit_trailingTangent;
        end
    end
    
    methods (Static)
        % TWO LINE INTERSECTIONS (FROM agent_VO)
        function [p_inter,isSuccessful] = findAny2DIntersection(P1,dP1,P2,dP2)
            % Find the intersection point between two 2D vectors. This
            % function isnt' interested if vertices are infront or behind
            % of the starting point.
            
            assert(numel(P1) == 2,'Input must be 2D');
            assert(numel(P2) == 2,'Input must be 2D');
            
            % SOME SUFFICIENTLY SMALL VALUE FOR INTERSECTION
            isSuccessful = logical(false);   % Default to no intersection
            p_inter = NaN(2,1);              % Default to no intersection
            
            % THE 2D DETERMININANT
            div = dP1(2)*dP2(1) - dP1(1)*dP2(2);
            if div == 0
                return % Lines are parallel
            end
            % THE SCALAR PROJECTIONS
            mua = (dP2(1)*(P2(2) - P1(2)) + dP2(2)*(P1(1) - P2(1))) / div;
%             mub = (dP1(1)*(P2(2) - P1(2)) + dP1(2)*(P1(1) - P2(1))) / div;
            % THE INTERSECTION POINT
            p_inter = P1 + mua*dP1;
            isSuccessful = logical(true);
        end
    end
end
% AGENT STATE VECTOR [x;y;phi;xdot;ydot;phidot]