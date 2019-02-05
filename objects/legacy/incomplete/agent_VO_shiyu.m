
classdef agent_VO_shiyu < agent_2D_VO
    
    %% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
    end
    %%  CLASS METHODS
    methods
        %% CONSTRUCTOR
        function obj = agent_VO_shiyu(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D_VO(varargin);                                        % Get super class 'agent'
            % DEFINE SOME KINEMATIC LIMITS
            obj.linearVelocityLimits = [2;2]; % Limit velocity in x/y
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        
        %% AGENT MAIN CYCLE
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});
            
            % /////////////////// WAYPOINT TRACKING ///////////////////////
            desiredSpeed = obj.linearVelocityLimits(1);
            desiredVelocity = [1;0]*desiredSpeed;
            if ~isempty(obj.targetWaypoint)
                waypointPosition = obj.targetWaypoint.state(1:2);
                desiredVelocity  = (waypointPosition/norm(waypointPosition))*desiredSpeed; % Desired relative velocity
            end
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            avoidanceSet = [agentSet,obstacleSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(desiredVelocity,avoidanceSet)
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            [d_speed,d_heading] = obj.simpleController(desiredVelocity);
            newHeadingRate = d_heading/dt;
            newSpeed = norm(obj.localState(7:8)) + d_speed;
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            newState = obj.stateDynamics_singleIntegrator(dt,newSpeed,newHeadingRate);
            
            % /// IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT /////////
%             if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
%                 newState = obj.freezeAgent();
%             end
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [newSpeed;newState(4:6)];         % Record the control inputs
        end

        % THE AVOIDANCE ROUTINE AT TIME (t=k)
        function [avoidanceVelocity] = getAvoidanceCorrection(obj,desiredVelocity,avoidanceSet)
            
            % OUR INFORMATION
            pa = obj.localState(1:2,1);            % Relative to ourselves
            va = obj.localState(7:8,1); % x,y velocity
            ra = obj.VIRTUAL.size;      % Our own radius
            
            VOset = [];
            for obstacle = 1:numel(avoidanceSet)
                pb = avoidanceSet(obstacle).state(1:2,1);
                vb = avoidanceSet(obstacle).state(3:4,1);
                rb = avoidanceSet(obstacle).size;
                
%                 [HRVO] = obj.defineHybridReciprocolVelocityObstacle(pa,va,ra,pb,vb,rb);
%                 [RVO]  = obj.defineReciprocolVelocityObstacle(pa,va,ra,pb,vb,rb);
                [VO]   = obj.defineVelocityObstacle(pa,va,ra,pb,vb,rb);
                VOset = horzcat(VOset,VO);
            end
            % COMPUTE THE OPTIMAL VELOCITY USING THE CLEAR PATH METHOD
            [optimalVelocity] = obj.clearPathMethod(va,desiredVelocity,VOset);
            
            avoidanceVelocity = optimalVelocity;
        end
        
        % CLEAR PATH METHOD
        function [optimalVelocity] = clearPathMethod(obj,va,va_prefer,VOset)
            % This function uses the clear path method to select the 
            % optimal velocity from the VO field.
            
            % GET THE NUMBER OF VELOCITY OBSTACLES
            VOsetLength = numel(VOset);
            if VOsetLength == 0
                optimalVelocity = va_prefer;
                return
            end
            % CALCULATE THE ORTHOGANOL PROJECTION POINTS
            dimensions = 2;
            projectionNum = 0;
            projectionPoints = zeros(dimensions,2*VOsetLength); % there will be 2*VONum projection points
            
            for VOnum = 1:numel(VOset)
                % MOVE THROUGH THE VELOCITY OBSTACLE SET
                VO_j = VOset(VOnum); 
                % left projection point
                [projPoint(:,1),isOnTheRay(1)] = obj.pointProjectionToRay(va,VO_j.apex,VO_j.edgeLeftUnit);
                % right projection point
                [projPoint(:,2),isOnTheRay(2)] = obj.pointProjectionToRay(va,VO_j.apex,VO_j.edgeRightUnit);
                % next check if every point is inside any VO; if yes, discard it 
                for projPointIdx = 1:dimensions
                    
                	if isOnTheRay(projPointIdx) == 1                       % If the projection point is on the ray
                        isPointOutsideAllVO = 1;                           % Initial value
                        
                        for idx = 1:VOsetLength
                            
                            if idx ~= VOnum
                                
                                % CHECK IF PROJECTION POINT IS INSIDE THE VO
                                if obj.isInsideVO(projPoint(:,projPointIdx),VO_j) == 1
                                    % if the point is inside any VO, then discard this point
                                    isPointOutsideAllVO = 0;
                                    break;
                                end
                            end
                            
                        end
                        
                        % CHECK IF ALL POINTS WERE ALL OUTSIDE THE VO
                        if isPointOutsideAllVO == 1
                            projectionNum = projectionNum + 1;
                            projectionPoints(:,projectionNum) = projPoint(:,projPointIdx);
%                           plot(projPoint(1,projPointIdx)+pa(1),projPoint(2,projPointIdx)+pa(2),'o','MarkerFaceColor','m');
                        end
                    end
                end   
            end
            
            % ///////////// CALCULATE THE INTERSECTION POINTS /////////////
            intersectionNum = 0;
            intersectionPoints = zeros(dimensions,4*VOsetLength*(VOsetLength-1)/2); % there will be at most 4*VONum*(VONum-1)/2 intersection points
            if VOsetLength >= 2 % there must be at least two VOs to calculate intersection points
                % ASSESS INTERSECTIONS FOR EACH VOi, VOj
                for VOIdx1 = 1:VOsetLength
                    VO_i = VOset(VOIdx1);
                    for VOIdx2 = (VOIdx1+1):VOsetLength
                        VO_j = VOset(VOIdx2);
                        % intersection between leftedge and leftedge
                        [p(:,1),isSuccessful(1)] = obj.twoRayIntersection(VO_i.apex,...
                                                                          VO_i.edgeLeftUnit,...
                                                                          VO_j.apex,...
                                                                          VO_j.edgeLeftUnit);
                        % intersection between leftedge and rightedge
                        [p(:,2),isSuccessful(2)] = obj.twoRayIntersection(VO_i.apex,...
                                                                          VO_i.edgeLeftUnit,...
                                                                          VO_j.apex,...
                                                                          VO_j.edgeRightUnit);
                        % intersection between rightedge and leftedge
                        [p(:,3),isSuccessful(3)] = obj.twoRayIntersection(VO_i.apex,...
                                                                          VO_i.edgeRightUnit,...
                                                                          VO_j.apex,...
                                                                          VO_j.edgeLeftUnit);
                        % intersection between rightedge and rightedge
                        [p(:,4),isSuccessful(4)] = obj.twoRayIntersection(VO_i.apex,...
                                                                          VO_i.edgeRightUnit,...
                                                                          VO_j.apex,...
                                                                          VO_j.edgeRightUnit);
                        
                        % next check if every point is inside any VO; if yes, discard it
                        % For each of VO's intersections
                        for intersectionIdx = 1:4
                            % if the intersection was valid
                            if isSuccessful(intersectionIdx) == 1
                                
                                isIntersectionPointOutsideVO = 1;          % initial value; will be changed later
                                
                                for VOIdx3 = 1:VOsetLength                 % ALL VO OBJECTS
                                    VO_k = VOset(VOIdx3);
                                    if VOIdx3 ~= VOIdx1 && VOIdx3 ~= VOIdx2
%                                         if obj.isInsideVO(p(:,intersectionIdx), VO_k)==1
                                        if obj.isInsideVO(p(:,intersectionIdx),VO_k) == 1
                                            isIntersectionPointOutsideVO = 0; % if the point is inside any VO, then discard this point
                                            break;
                                        end
                                    end
                                    
                                end
                                if isIntersectionPointOutsideVO==1
                                    intersectionNum=intersectionNum+1;
                                    intersectionPoints(:,intersectionNum)=p(:,intersectionIdx);
                                    %                                     plot(intersectionPoints(1,intersectionNum)+pa(1),intersectionPoints(2,intersectionNum)+pa(2),'o','MarkerFaceColor','g');
                                end
                            end
                        end
                    end
                end
                
            end
            % ///// CHOOSE OPTIMAL VELOCITY FROM THE CANDIDATE POINTS /////
            optimalityMetricDistance = inf;
            optimalVelocity = zeros(2,1); % if no feasible, then the value is unchanged as zero
            compareVelocity=va_prefer;
            if projectionNum + intersectionNum>0
                candidatePoints = [projectionPoints(:,1:projectionNum),intersectionPoints(:,1:intersectionNum)];
                for k = 1:size(candidatePoints,2)
                    dis = norm(compareVelocity-candidatePoints(:,k));
                    if dis<optimalityMetricDistance
                        optimalVelocity = candidatePoints(:,k);
                        optimalityMetricDistance = dis;
                    end
                end
            else
                disp('there is no feasible velocity!')
            end
        end
        % ORCA NEW VELOCITY METHOD
        function [optimalVelocity] = newVelocityORCAMethod(VOset, va, va_prefer)
            
            dimensions = 2;
            VOnum = numel(VOset);
            
            % calcualte A and b for quadratic programing: A*x>=b
            quadraticA = zeros(VOnum,dimensions);
            quadraticb = zeros(VOnum,1);
            
            % FOR EACH VELOCITY OBSTACLE
            for k = 1:VOnum
                % COLLECT THE PROJECTION POINTS FOR EACH OBSTACLE
                [projPoint(:,1), isOnTheRay(1)] = obj.pointProjectionToRay(va,VOset(k).apex,VOset(k).edgeRightUnit);
                [projPoint(:,2), isOnTheRay(2)] = obj.pointProjectionToRay(va,VOset(k).apex,VOset(k).edgeLeftUnit);
                if isOnTheRay(1) == 1 && isOnTheRay(2) == 1
                    if norm(projPoint(:,1)-va) <= norm(projPoint(:,2)-va)
                        % va is close to the right edge
                        normalVectorHalfPlane=projPoint(:,1) - va;
                        pointHalfPlane=(va+projPoint(:,1))/2;              % agent A takes half of the responsibility;
                    else
                        % va is close to the left edge
                        normalVectorHalfPlane = projPoint(:,2) - va;
                        pointHalfPlane=(va+projPoint(:,2))/2;              % agent A takes half of the responsibility;
                    end
                    
                elseif isOnTheRay(1) == 1
                    normalVectorHalfPlane=projPoint(:,1)-va;
                    pointHalfPlane=(va+projPoint(:,1))/2;                  % agent A takes half of the responsibility;
                else
                    normalVectorHalfPlane=projPoint(:,2)-va;
                    pointHalfPlane=(va+projPoint(:,2))/2;                  % agent A takes half of the responsibility;
                end
                quadraticA(k,:) = normalVectorHalfPlane';                  %/norm(normalVectorHalfPlane);
                quadraticb(k)=normalVectorHalfPlane'*pointHalfPlane;
            end
            % COMPUTE THE OPTIMAL VELOCITY FROM THE SOLUTION SET
            [optimalVelocity, fval, exitFlag] = quadprog(2*eye(dimensions),...
                                               -2*va_prefer,...
                                               -quadraticA,...
                                               -quadraticb);
        end
        
        % DEFINE TRUNCATED VELOCITY OBSTACLE
        function [trunVO] = defineTruncatedVelocityObstacle(pa,va,ra,pb,vb,rb,tau)
            % This function builds the truncated velocity obstacles.
            
            % GET THE ORIGINAL VELOCITY OBSTACLE
            trunVO = obj.defineVelocityObstacle(pa,va,ra,pb,vb,rb);
            % TRUNCATE THE VO IF 'va' IS INSIDE VO
            if trunVO.isVaInsideCone
               trunVO.truncationTau=tau;
               trunVO.truncationCircleCenter=(pb-pa)/tau + vb;
               trunVO.truncationCircleRadius=(ra+rb)/tau;
            end
        end
        % DEFINE THE HYBRID RECIPROCOL VELOCITY OBSTACLE
        function [HRVO] = defineHybridReciprocolVelocityObstacle(obj,pa,va,ra,pb,vb,rb)
            % This function generates the hybrid reciprocol velocity 
            % obstacle by calculating the new apex position from the VO and
            % RVO geometry.
            
            % GET THE VELOCITY OBSTACLE
            VO  = obj.defineVelocityObstacle(pa,va,ra,pb,vb,rb);
            % GET THE RECIPROCOL VELOCITY OBSTACLE
            RVO = obj.defineReciprocolVelocityObstacle(pa,va,ra,pb,vb,rb);
            
            % COMPARE THE ORIENTATION OF THE VO
            if VO.isVaRightCenterLine==1
                p1 = VO.apex;           % VO apex
                v1 = VO.edgeLeftUnit;   % VO edge
                p2 = RVO.apex;          % RVO apex
                v2 = RVO.edgeRightUnit; % RVO edge
            else
                p1 = VO.apex;           % VO apex
                v1 = VO.edgeRightUnit;  % VO edge
                p2 = RVO.apex;          % RVO apex
                v2 = RVO.edgeLeftUnit;  % RVO edge
            end
            
            if norm(v2-v1*v1'*v2/(v1'*v1)) == 0 || norm(v2-v1*v1'*v2/(v1'*v1)) < 0.0001
                % v1 and v2 are parallel: this happens when norm(pa-pb)=ra+rb; v1 and
                % v2 are parallel and pointing to reverse directions
                % in this case, simply choose HRVO as RVO
                HRVOapex = RVO.apex;
                %     disp('warning calculate HRVO!')
            else
                Proj_v1 = (eye(2)-v1*v1'/(v1'*v1));
                Proj_v2 = (eye(2)-v2*v2'/(v2'*v2));
                HRVOapex = (Proj_v1+Proj_v2)\(Proj_v1*p1+Proj_v2*p2); % here \ means inv(A)*b
            end
            HRVO = RVO;
            HRVO.apex = HRVOapex;
        end
        % DEFINE THE RECIPRICOL VELOCITY OBSTACLE
        function [RVO]  = defineReciprocolVelocityObstacle(obj,pa,va,ra,pb,vb,rb)
            % This function generates the reciprocol velocity obstacle from
            % the original velocity obstacle.
            
            % GET THE VELOCITY OBSTACLES
            [RVO] = obj.defineVelocityObstacle(pa,va,ra,pb,vb,rb);
            % AUGMENT THE APEX OF THE VO
            RVO.apex = (va + vb)/2;
        end
        % DEFINE THE VELOCITY OBSTACLE
        function [VO] = defineVelocityObstacle(obj,pa,va,ra,pb,vb,rb)
            % RELATIVE POSITION PARAMETERS
            centerLine = pb - pa;                   % Relative position
            norm_centerLine = norm(centerLine);
            unit_centerLine = centerLine/norm_centerLine;
            
            % GET THE OPEN ANGLE
            if norm_centerLine <= (ra+rb) 
                % special case: alreayd collided!
                halfOpenAngle=pi/2;
            else
                halfOpenAngle=asin((ra+rb)/norm_centerLine);
            end
            openAngle=2*halfOpenAngle;
            
            % GET THE PLANAR VERTICES 
            % These were moved out to allow the generation of HRVO using
            % the same function.
            unit_edgeLeft = obj.rotateCounterclockwise(unit_centerLine,halfOpenAngle);
            unit_edgeRight = obj.rotateCounterclockwise(unit_centerLine,-halfOpenAngle);
            
            % ASSEMBLE INITIAL VO PARAMETERS
            VO = struct('apex',vb,...
                   'openAngle',openAngle,...
              'centerLineUnit',unit_centerLine,...                
                'edgeLeftUnit',unit_edgeLeft,...
               'edgeRightUnit',unit_edgeRight); 
              
          
            % DETERMINE IF 'va' IS INSIDE THE VO
            if obj.isInsideVO(va,VO)
                isVaInsideCone = 1;
                if unit_edgeLeft'*(va-vb) < unit_edgeRight'*(va-vb) % the angle to the left edge is large then the angle to the right edge
                    isVaRightCenterLine=1;
                else
                    isVaRightCenterLine=0;
                end     
            else
                isVaInsideCone = 0;
                isVaRightCenterLine = NaN; % If va is outside the VO, other parameters are NaN.
            end
            % ADDITIONAL PLACEHOLDER VO PARAMETERS
            VO.isVaInsideCone = logical(isVaInsideCone);
            VO.isVaRightCenterLine = isVaRightCenterLine; % 0,1,NaN
            VO.truncationTau = 0;
            VO.truncationCircleCenter = [0;0];
            VO.truncationCircleRadius = 0;
        end
    end
    % STATIC VELOCITY OBSTACLE METHODS
    methods (Static)        
        % CHECK THE AGENT VELOCITY IS INSIDE THE VO      
        function [flag] = isInsideVO(point, VO)
        % determine if the point p is inside the given VO
        % angle<halfOpenAngle; here numericalTolerance is an error tolarance
        numericalTolerance=-(1e-8);
        flag = 0;
        if VO.centerLineUnit'*(point-VO.apex) - norm(point-VO.apex)*cos(VO.openAngle/2) > numericalTolerance
            flag=1;
        end
        end
        
        % GET THE TWO RAY INTERSECTION
        function [p, isSuccessful] = twoRayIntersection(p1,v1,p2,v2)
        % p1 is one point on line 1, v1 is a vector parallel to line 1
        % first of all, if v1 and v2 are parallel, then there is either no
        % intersection or infinite intersection points. In this case do not
        % calculate and directly return
        
        if norm(v2-v1*v1'*v2/(v1'*v1))==0 % v1 and v2 are parallel
            isSuccessful=0;
            p = [NaN, NaN]';
        else
            Proj_v1 = (eye(2)-v1*v1'/(v1'*v1));
            Proj_v2 = (eye(2)-v2*v2'/(v2'*v2));
            p=(Proj_v1+Proj_v2)\(Proj_v1*p1+Proj_v2*p2); % here \ means inv(A)*b
            % the intersection point must be on the ray (not the other half ray)
            if v1'*(p-p1)>0 && v2'*(p-p2)>0 % on the ray
                isSuccessful = 1;
            else
                isSuccessful = 0;
                p=[NaN, NaN]';
            end
        end
        end
        % GET THE POINT PROJECTION ON RAY
        function [ projectedPoint,isOnTheRay] = pointProjectionToRay(p,p0,v0)
            % INPUTS:
            % p  - Is the point to be projected.
            % p0,v0 - The line defining points
            % OUTPUTS:
            % 
            projectedPoint=v0*v0'/(v0'*v0)*(p-p0)+p0;

            if v0'*(projectedPoint-p0)>0 % if on the ray
                isOnTheRay = 1;
            else
                isOnTheRay = 0;
            end
        end
        % ROTATE 2D VECTOR COUNTER-CLOCKWISE
        function vNew = rotateCounterclockwise(vOld,angle)
            R = [cos(angle),-sin(angle);
                 sin(angle),cos(angle)];
            vNew = R*vOld;
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]