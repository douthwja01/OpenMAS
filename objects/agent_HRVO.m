%% HYBRID RECIPROCAL VELOCITY OBSTACLE METHOD (agent_HRVO.m) %%%%%%%%%%%%%%
% This agent classe implements the Hybrid-Reciprocal Velocity Obstacle
% (HRVO) method outlined in: "The Hybrid Reciprocal Velocity Obstacle" by
% Jur van den Berg et al.

% Author: James A. Douthwaite

classdef agent_HRVO < agent_RVO
    properties
        % PROPERTIES UNIQUE TO THE HRVO METHOD
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function obj = agent_HRVO(varargin)
            
            % Call the super class
            obj@agent_RVO(varargin);  
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % ///////////////////////////////////////////////////////////// 
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.

        % GET THE VO VELOCITY CORRECTION
        function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,visualiseProblem)
            % This function calculates the collision avoidance heading and
            % velocity correction.
            
            % AGENT KNOWLEDGE
            [p_i,v_i,r_i] = obj.GetAgentMeasurements();                    % Its own position, velocity and radius
            % GET OBSTACLE DATA
            obstacleIDs  = [obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.obstacle).objectID];
            agentIDs     = [obj.MEMORY([obj.MEMORY.type] == OMAS_objectType.agent).objectID];
            avoidanceIDs = [agentIDs,obstacleIDs];
            
            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:numel(avoidanceIDs)
                % Obstacle data
                p_j = obj.GetLastMeasurementByObjectID(avoidanceIDs(item),'position');
                
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;            % Maximum number of neighbours
                neighbourConditionB = norm(p_j) < obj.neighbourDist;       % [CONFIRMED] 
                if ~neighbourConditionA || ~neighbourConditionB
                    continue
                end           
                % More obstacle information
                v_j = obj.GetLastMeasurementByObjectID(avoidanceIDs(item),'velocity');
                r_j = obj.GetLastMeasurementByObjectID(avoidanceIDs(item),'radius');                                
                % Make none-relative
                p_j = p_j + p_i; 
                v_j = v_j + v_i;                                           % Convert relative parameters to absolute
                tau_j = 0;

                % OBSTACLE TYPE BEHAVIOUR
                if OMAS_objectType.agent == obj.GetLastMeasurementByObjectID(avoidanceIDs(item),'type')
                    % DEFINE RVO AS AGENT BEHAVIOUR
                    [VO_i] = obj.define3DHybridReciprocalVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                else
                    % OBSTACLE BEHAVIOUR
                    [VO_i] = obj.define3DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                end       
                % CONCATINATE THE VO SET
                VO = [VO,VO_i]; 
            end
            
            % GET THE CAPABLE VELOCITIES
            capableVelocities = obj.feasabliltyMatrix;
            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
            escapeVelocities  = obj.getEscapeVelocities(capableVelocities,VO);              % The viable velocities from the capable 
            % APPLY THE MINIMUM DIFFERENCE SEARCH STRATEGY
            avoidanceVelocity = obj.strategy_minimumDifference(desiredVelocity,escapeVelocities);                  
            
            % SPECIAL CASE- VELOCITY MAGNITUDE IS ZERO
            speed = norm(avoidanceVelocity);
            headingVector = avoidanceVelocity/speed;
            if isnan(headingVector)
                headingVector = [1;0;0];  % Retain previous heading
            end
            
            % PLOT THE VECTOR CONSTRUCT
            if visualiseProblem && obj.objectID == visualiseProblem
                % PLOT THE LEADING TANGENT VECTOR
                OMAS_axisTools.drawTriad(p_i,eye(3));
                
                % CURRENT VELOCITY
                q = quiver3(gca,p_i(1),p_i(2),p_i(3),v_i(1),v_i(2),v_i(3),'m');
                q.AutoScaleFactor = 1;
                % DESIRED VELOCITY
                q = quiver3(gca,p_i(1),p_i(2),p_i(3),desiredVelocity(1),desiredVelocity(2),desiredVelocity(3),'g');
                q.AutoScaleFactor = 1;   
                % AVOIDANCE VELOCITY
                q = quiver3(gca,p_i(1),p_i(2),p_i(3),avoidanceVelocity(1),avoidanceVelocity(2),avoidanceVelocity(3),'b');
                q.AutoScaleFactor = 1;
                
%                 plotCapableVelocities = capableVelocities + p_i;
%                 plotEscapeVelocities = escapeVelocities + p_i;
%                 plotAvoidanceVelocities = avoidanceVelocity + p_i;
%                 scatter3(gca,plotCapableVelocities(1,:),plotCapableVelocities(2,:),plotCapableVelocities(3,:),'r'); 
%                 scatter3(gca,plotEscapeVelocities(1,:),plotEscapeVelocities(2,:),plotEscapeVelocities(3,:),'g');                 
%                 scatter3(gca,plotAvoidanceVelocities(1),plotAvoidanceVelocities(2),plotAvoidanceVelocities(3),'b','filled'); 
            end
        end
    end
    
    % /////////// HYBRID RECIPRICOL VELOCITY OBSTACLE METHODS /////////////
    methods (Access = public)
        % DEFINE THE 3D HYBRID RECIPCOCAL VELOCITY OBSTACLE (HRVO)
        function [HRVO] = define3DHybridReciprocalVelocityObstacle(obj,p_i,v_i,r_i,p_j,v_j,r_b,tau_b,plotOn)
            % This function designs a Hybrid-Reciprocal Velocity Obstacle
            % (HRVO) for a sensed obstacle. This method is that presented
            % in the paper "The Hybrid Reciprocal Velocity Obstacle" by
            % Jamie Snape et. al.
            % INPUTS:
            % obj      - The agent object
            % obstacle - The sensed obstacle from the agents memory
            % plotOn   - The output plot flag
            % OUTPUTS:
            % VOapex   - The apex of HRVO
            % VOaxis   - The vector describing the cone axis
            % gamma_VO - The cone angle from the axis vector

            % GENERATE A STANDARD VELOCITY OBSTACLE
            [VO]  = obj.define3DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_b,tau_b,plotOn);
            [RVO] = obj.define3DReciprocalVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_b,tau_b,plotOn);       
            
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
%             [HRVOapex, isSuccessful] = obj.twoRayIntersection2D(VO.apex,VOelement,RVO.apex,-RVOelement); % Project the RVO vector towards the VO element
%             [pa,pb,isSuccessful] = obj.twoRayIntersection3D(P1,dP1,P3,dP3)
%             [HRVOapex,pb,isSuccessful] = obj.findAny3DIntersection(VO.apex,VOelement,RVO.apex,-RVOelement);
            [HRVOapex,pb,isSuccessful] = GetRayRayIntersection3D_mex(VO.apex,VOelement,RVO.apex,-RVOelement);
            
            if any(isnan(HRVOapex))%~isSuccessful
                HRVOapex = RVO.apex;
                warning('no line intersection');
            end
                        
            % //////// CONSTRUCT THE HRVO FROM THE RVO TEMPLATE ///////////
            HRVO = RVO;
            HRVO.apex = HRVOapex;
            HRVO.leadingEdgeUnit = unit_leadingTangent;
            HRVO.trailingEdgeUnit = unit_trailingTangent;
            
            % PROBLEM VISUALISATION
            if plotOn && obj.objectID == plotOn
                % DRAW THE VO
                coneApex   = p_i + HRVO.apex;
                coneCenter = coneApex + HRVO.axisUnit*HRVO.axisLength;
                % DRAW OBSTACLE AS INPUTTED
                [Xb,Yb,Zb] = obj.getSphere(p_j,r_b);
                sphB = mesh(gca,Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                    'FaceColor','b',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                % OBSTACLE VELOCITY (ABS)
                q = quiver3(gca,p_j(1),p_j(2),p_j(3),...
                    v_j(1),v_j(2),v_j(3),'b','filled');          % Plot the obstacles velocity from its relative position
                q.AutoScaleFactor = 1;
                
                % DRAW THE VELOCITY OBSTACLE PROJECTION
                [X_vo,Y_vo,Z_vo] = obj.getSphere(coneCenter,(r_i + r_b));
                sphB = mesh(gca,X_vo,Y_vo,Z_vo);
                set(sphB,'facealpha',0.2,...
                    'FaceColor','b',...
                    'LineWidth',0.1,...
                    'EdgeAlpha',0);
                
                % //////////////// VO CONSTRUCTION ////////////////////////
                % THE LEADING TANGENT                
                leadingTangent = HRVO.leadingEdgeUnit*HRVO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                                leadingTangent(1),leadingTangent(2),leadingTangent(3),'b');                  
                q.AutoScaleFactor = 1;
                
                % THE TRAILING TANGENT                
                trailingTangent = HRVO.trailingEdgeUnit*HRVO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                            trailingTangent(1),trailingTangent(2),trailingTangent(3),'b');
                q.AutoScaleFactor = 1;
                
                % ASSEMBLE VELOCITY OBSTACLE CONES
                try
                    tangentPoint = coneApex + leadingTangentVector;
                    nodes = 11;
                    obj.vectorCone(coneApex,coneCenter,tangentPoint,nodes,'b');
                catch
                end
            end
        end
    end
    
    methods (Static)
        % GET THE INTERSECTION BETWEEN TWO LINES IN 3D
        function [pa,pb,isSuccessful] = findAny3DIntersection(P1,dP1,P3,dP3)
            %     /*
            %        Calculate the line segment PaPb that is the shortest route between
            %        two lines P1P2 and P3P4. Calculate also the values of mua and mub where
            %           Pa = P1 + mua (P2 - P1)
            %           Pb = P3 + mub (P4 - P3)
            %        Return FALSE if no solution exists.
            %     */
            
            % SOME SUFFICIENTLY SMALL VALUE FOR INTERSECTION
            EPS = 1E-9;
            isSuccessful = logical(false);
            pa = NaN(3,1);
            pb = NaN(3,1);
            
            % GET THE DIRECTIONS
            p13(1) = P1(1) - P3(1);
            p13(2) = P1(2) - P3(2);
            p13(3) = P1(3) - P3(3);
            
%             if (abs(dP3(1)) < EPS && abs(dP3(2)) < EPS && abs(dP3(3)) < EPS)
%                 return
%             end
%             
%             if (abs(dP1(1)) < EPS && abs(dP1(2)) < EPS && abs(dP1(3)) < EPS)
%                 return
%             end
            
            d1343 = p13(1) * dP3(1) + p13(2) * dP3(2) + p13(3) * dP3(3);
            d4321 = dP3(1) * dP1(1) + dP3(2) * dP1(2) + dP3(3) * dP1(3);
            d1321 = p13(1) * dP1(1) + p13(2) * dP1(2) + p13(3) * dP1(3);
            d4343 = dP3(1) * dP3(1) + dP3(2) * dP3(2) + dP3(3) * dP3(3);
            d2121 = dP1(1) * dP1(1) + dP1(2) * dP1(2) + dP1(3) * dP1(3);
            
            denom = d2121 * d4343 - d4321 * d4321;
            if (abs(denom) < EPS)
                return
            end
            
            numer = d1343 * d4321 - d1321 * d4343;
            
            mua = numer / denom;
            mub = (d1343 + d4321 * mua) / d4343;
            
            % BUILD INTERSECTION COORDINATES
            pa(1) = P1(1) + mua * dP1(1);
            pa(2) = P1(2) + mua * dP1(2);
            pa(3) = P1(3) + mua * dP1(3);
            pb(1) = P3(1) + mub * dP3(1);
            pb(2) = P3(2) + mub * dP3(2);
            pb(3) = P3(3) + mub * dP3(3);
            % IS SUCCESSFUL
            isSuccessful = logical(true);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]