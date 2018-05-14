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
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_RVO(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent_VO(varargin);                                        % Get the supercalss

            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        
        % MAIN CYCLE IS INHERITED FROM THE SUPERCLASS
        % The waypoint handling, dynamics and update procedures must also
        % be the same as the super-class.    
        
        % GET THE VO VELOCITY CORRECTION
        function [headingVector,speed] = getAvoidanceCorrection(obj,desiredVelocity,knownObstacles,visualiseProblem)
            % This function calculates the collision avoidance heading and
            % velocity correction.
            
            % AGENT KNOWLEDGE
            p_a = obj.localState(1:3,1);            % In relative frame centered at a
            v_a = obj.localState(7:9,1);            % Body axis velocity
            r_a = obj.VIRTUAL.radius;

            % BUILD THE VELOCITY OBSTACLE SET
            VO = [];
            for item = 1:length(knownObstacles)
                % NEIGHBOUR CONDITIONS
                neighbourConditionA = item < obj.maxNeighbours;                                 % Maximum number of neighbours
                neighbourConditionB = norm(knownObstacles(item).position) < obj.neighbourDist;  % [CONFIRMED] 
                neighbourConditionC = ~any(isnan(knownObstacles(item).velocity));               % Wait for a valid velocity reading
                if ~neighbourConditionA || ~neighbourConditionB || ~neighbourConditionC
                    continue
                end
                % OBSTACLE KNOWLEDGE
                p_b = knownObstacles(item).position + p_a; 
                v_b = knownObstacles(item).velocity + v_a;                 % Convert relative parameters to absolute
                r_b = knownObstacles(item).radius;
                tau = 0;
                % DEFINE THE VELOCITY OBSTACLE PROPERTIES
                [RVO_i] = obj.define3DReciprocalVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau,visualiseProblem);
                VO = vertcat(VO,RVO_i);
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
                OMAS_axisTools.drawTriad(p_a,eye(3));
                
                % CURRENT VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),v_a(1),v_a(2),v_a(3),'m');
                q.AutoScaleFactor = 1;
                % DESIRED VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),desiredVelocity(1),desiredVelocity(2),desiredVelocity(3),'g');
                q.AutoScaleFactor = 1;   
                % AVOIDANCE VELOCITY
                q = quiver3(gca,p_a(1),p_a(2),p_a(3),avoidanceVelocity(1),avoidanceVelocity(2),avoidanceVelocity(3),'b');
                q.AutoScaleFactor = 1;
                
                plotCapableVelocities = capableVelocities + p_a;
                plotEscapeVelocities = escapeVelocities + p_a;
                plotAvoidanceVelocities = avoidanceVelocity + p_a;
                scatter3(gca,plotCapableVelocities(1,:),plotCapableVelocities(2,:),plotCapableVelocities(3,:),'r'); 
                scatter3(gca,plotEscapeVelocities(1,:),plotEscapeVelocities(2,:),plotEscapeVelocities(3,:),'g');                 
                scatter3(gca,plotAvoidanceVelocities(1),plotAvoidanceVelocities(2),plotAvoidanceVelocities(3),'b','filled'); 
            end
        end
    end
    
    % /////////////// RECIPRICOL VELOCITY OBSTACLE METHODS ////////////////
    methods (Access = public)
        % DEFINE THE 3D RECIPROCAL VELOCITY OBSTACLE (RVO)
        function [RVO] = define3DReciprocalVelocityObstacle(obj,p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn)
            % This function defines a reciprocol velocity obstacle using
            % the available obstacle knowledge.
            
            
            % GENERATE A STANDARD VELOCITY OBSTACLE
            [VO] = obj.define3DVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau_b,plotOn);
             
            % ////// DEFINE RECIPROCAL VELOCITY OBSTACLE PARAMETERS ///////
            RVO = struct('apex',(v_a + v_b)/2,...
                     'axisUnit',VO.axisUnit,...
                   'axisLength',VO.axisLength,...
                    'openAngle',VO.openAngle,...  
              'leadingEdgeUnit',VO.leadingEdgeUnit,...
             'trailingEdgeUnit',VO.trailingEdgeUnit,...
                  'isVaLeading',VO.isVaLeading); 
            % ///////////////////////////////////////////////////////////// 
            % PROBLEM VISUALISATION
            if plotOn && obj.objectID == plotOn
                % DRAW THE VO
                coneApex   = p_a + RVO.apex;
                coneCenter = coneApex + RVO.axisUnit*RVO.axisLength; 
                % DRAW OBSTACLE AS INPUTTED
                [Xb,Yb,Zb] = obj.getSphere(p_b,r_b);                            
                sphB = mesh(gca,Xb,Yb,Zb);  % Obstacle
                set(sphB,'facealpha',0.2,...
                         'FaceColor','b',...
                         'LineWidth',0.1,...
                         'EdgeAlpha',0); 
                % OBSTACLE VELOCITY (ABS)
                q = quiver3(gca,p_b(1),p_b(2),p_b(3),...
                                v_b(1),v_b(2),v_b(3),'b','filled');          % Plot the obstacles velocity from its relative position
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