%% RECIPRICOL VELOCITY OBSTACLE METHOD (agent_RVO.m) %%%%%%%%%%%%%%%%%%%%%%%
% This agent implements the reciprocal velocity obstacle (RVO) algorithm
% implemented in " Reciprocal Velocity Obstacles for real-time multi-agent 
% navigation - Jur van den Berg et al 2008"

% Author: James A. Douthwaite

classdef agent_RVO < agent_VO
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTIES UNIQUE TO THE RVO METHOD
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function obj = agent_RVO(varargin)
            % Get the superclass
            obj@agent_VO(varargin);  

            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.    
        
        % GET THE VO VELOCITY CORRECTION
        function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,knownObstacles,visualiseProblem)
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
                    [VO_i] = obj.define3DReciprocalVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                else
                    % OBSTACLE BEHAVIOUR
                    [VO_i] = obj.define3DVelocityObstacle(p_i,v_i,r_i,p_j,v_j,r_j,tau_j,visualiseProblem);
                end       
                % CONCATINATE THE VO SET
                VO = [VO,VO_i]; 
            end
            
            % GET THE CAPABLE VELOCITIES
            capableVelocities = obj.feasabliltyMatrix;                     % Get all the feasible velocities this timestep 
            
            % SUBTRACT THE VELOCITY OBSTACLES FROM THE VELOCITY FIELDS
            [escapeVelocities] = obj.getEscapeVelocities(capableVelocities,VO);              % The viable velocities from the capable 
            
            % APPLY THE MINIMUM DIFFERENCE SEARCH STRATEGY
            [avoidanceVelocity] = obj.strategy_minimumDifference(desiredVelocity,escapeVelocities);                  
            
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
                
                plotCapableVelocities = capableVelocities + p_i;
                plotEscapeVelocities = escapeVelocities + p_i;
                plotAvoidanceVelocities = avoidanceVelocity + p_i;
                scatter3(gca,plotCapableVelocities(1,:),plotCapableVelocities(2,:),plotCapableVelocities(3,:),'r'); 
                scatter3(gca,plotEscapeVelocities(1,:),plotEscapeVelocities(2,:),plotEscapeVelocities(3,:),'g');                 
                scatter3(gca,plotAvoidanceVelocities(1),plotAvoidanceVelocities(2),plotAvoidanceVelocities(3),'b','filled'); 
            end
        end
    end
    
    % /////////////// RECIPRICOL VELOCITY OBSTACLE METHODS ////////////////
    methods (Access = public)
        % DEFINE THE 3D RECIPROCAL VELOCITY OBSTACLE (RVO)
        function [RVO] = define3DReciprocalVelocityObstacle(obj,p_i,v_i,r_a,p_j,v_j,r_b,tau_b,plotOn)
            % This function defines a reciprocol velocity obstacle using
            % the available obstacle knowledge.
            
            
            % GENERATE A STANDARD VELOCITY OBSTACLE
            [RVO] = obj.define3DVelocityObstacle(p_i,v_i,r_a,p_j,v_j,r_b,tau_b,plotOn);
            % ////// DEFINE RECIPROCAL VELOCITY OBSTACLE PARAMETERS ///////
            RVO.apex = (v_i + v_j)/2; % Modify the apex position

            % ///////////////////////////////////////////////////////////// 
            % PROBLEM VISUALISATION
            if plotOn && obj.objectID == plotOn
                % DRAW THE VO
                coneApex   = p_i + RVO.apex;
                coneCenter = coneApex + RVO.axisUnit*RVO.axisLength; 
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
                [X_vo,Y_vo,Z_vo] = obj.getSphere(coneCenter,(r_a + r_b));  
                sphB = mesh(gca,X_vo,Y_vo,Z_vo); 
                set(sphB,'facealpha',0.2,...
                         'FaceColor','r',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0);
                
                % //////////////// VO CONSTRUCTION ////////////////////////
                % THE LEADING TANGENT                
                leadingTangent = VO.leadingEdgeUnit*VO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                                leadingTangent(1),leadingTangent(2),leadingTangent(3),'r');                  
                q.AutoScaleFactor = 1;
                
                % THE TRAILING TANGENT                
                trailingTangent = VO.trailingEdgeUnit*VO.axisLength;
                q = quiver3(gca,coneApex(1),coneApex(2),coneApex(3),...
                            trailingTangent(1),trailingTangent(2),trailingTangent(3),'r');
                q.AutoScaleFactor = 1;
                
                % ASSEMBLE VELOCITY OBSTACLE CONES
                try
                    tangentPoint = coneApex + leadingTangentVector;
                    nodes = 11;
                    obj.vectorCone(coneApex,coneCenter,tangentPoint,nodes,'r');
                catch
                end
            end
        end
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]