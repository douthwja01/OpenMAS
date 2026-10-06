%% VELOCITY OBSTACLE TOOLS (VO_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared VO geometry used by 2D and 3D avoidance agents.

classdef VO_tools
    methods (Static)
        function [VO] = defineVelocityObstacle(p_a,v_a,r_a,p_b,v_b,r_b,tau,epsilon)
            if nargin < 8 || isempty(epsilon)
                epsilon = 1E-5;
            end

            lambda_ab = p_b - p_a;
            r_c = r_a + r_b;
            mod_lambda_ab = norm(lambda_ab);
            unit_lambda_ab = lambda_ab/mod_lambda_ab;

            referenceAxis = [1;0;0];
            planarNormal = cross(referenceAxis,unit_lambda_ab);
            if sum(planarNormal) == 0
                planarNormal = cross([0;0;1],unit_lambda_ab);
            end

            halfAlpha = real(asin(r_c/mod_lambda_ab));
            [leadingTangentVector] = OMAS_geometry.rodriguesRotation(lambda_ab,planarNormal,halfAlpha);
            VOaxis = (dot(leadingTangentVector,lambda_ab)/mod_lambda_ab^2)*lambda_ab;
            axisLength = norm(VOaxis);
            axisUnit = VOaxis/axisLength;

            [leadingTangentVector] = OMAS_geometry.rodriguesRotation(VOaxis,planarNormal,halfAlpha);
            unit_leadingTangent = leadingTangentVector/norm(leadingTangentVector);
            [trailingTangentVector] = OMAS_geometry.rodriguesRotation(VOaxis,planarNormal,-halfAlpha);
            unit_trailingTangent = trailingTangentVector/norm(trailingTangentVector);

            isVaLeading = unit_leadingTangent'*(v_a - v_b) > unit_trailingTangent'*(v_a - v_b);

            VO = struct('apex',v_b,...
                'axisUnit',axisUnit,...
                'axisLength',mod_lambda_ab,...
                'openAngle',2*halfAlpha,...
                'leadingEdgeUnit',unit_leadingTangent,...
                'trailingEdgeUnit',unit_trailingTangent,...
                'isVaLeading',isVaLeading,...
                'isVaInsideCone',0,...
                'truncationTau',tau,...
                'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
                'truncationCircleRadius',(r_a + r_b)/tau);
            VO.isVaInsideCone = VO_tools.isInCone(VO,v_a,epsilon);
        end

        function [RVO] = reciprocalObstacle(VO,v_i,v_j)
            RVO = VO;
            RVO.apex = (v_i + v_j)/2;
        end

        function [VO] = projectTo2D(VO)
            fields = {'apex','axisUnit','leadingEdgeUnit','trailingEdgeUnit','truncationCircleCenter'};
            for i = 1:numel(fields)
                value = VO.(fields{i});
                if numel(value) > 2
                    VO.(fields{i}) = value(1:2,1);
                end
            end
        end

        function [flag] = isInCone(VO,probePoint,epsilon)
            if nargin < 3 || isempty(epsilon)
                epsilon = 1E-5;
            end
            probeVector = probePoint - VO.apex;
            unit_probeVector = probeVector/norm(probeVector);
            probeDot = dot(unit_probeVector,VO.axisUnit);
            theta = real(acos(probeDot));
            flag = 0;
            if (theta - VO.openAngle/2 < epsilon) && probeDot > 0
                flag = 1;
            end
        end

        function [VO] = defineComplexVelocityObstacle(p_a,v_a,r_a,p_b,v_b,geometry,tau)
            norm_lateralPositionProjection = norm([p_b(1:2);0]);
            unit_lateralPositionProjection = [p_b(1:2);0]/norm_lateralPositionProjection;

            a_min = 0;
            a_max = 0;
            unit_leadingTangent = unit_lateralPositionProjection;
            unit_trailingTangent = unit_lateralPositionProjection;
            for v = 1:size(geometry.vertices,1)
                lateralVertexProjection = [geometry.vertices(v,1:2)';0];
                unit_lateralVertexProjection = lateralVertexProjection/norm(lateralVertexProjection);
                crossProduct = cross(unit_lateralVertexProjection,unit_lateralPositionProjection);
                dotProduct = dot(unit_lateralVertexProjection,unit_lateralPositionProjection);
                signedAngularProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
                if signedAngularProjection > a_max
                    a_max = signedAngularProjection;
                    unit_trailingTangent = unit_lateralVertexProjection;
                elseif signedAngularProjection < a_min
                    a_min = signedAngularProjection;
                    unit_leadingTangent = unit_lateralVertexProjection;
                end
            end

            crossProduct = cross(v_a,unit_lateralPositionProjection);
            dotProduct = dot(v_a/norm(v_a),unit_lateralPositionProjection);
            signedVelocityProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
            isVaInsideCone = signedVelocityProjection < a_max && signedVelocityProjection > a_min;
            equivalentOpenAngle = abs(a_min) + abs(a_max);
            effectiveRadii = norm_lateralPositionProjection*sin(equivalentOpenAngle);

            VO = struct('apex',v_b,...
                'axisUnit',unit_lateralPositionProjection,...
                'axisLength',norm_lateralPositionProjection,...
                'openAngle',equivalentOpenAngle,...
                'leadingEdgeUnit',unit_leadingTangent,...
                'trailingEdgeUnit',unit_trailingTangent,...
                'isVaLeading',1,...
                'isVaInsideCone',isVaInsideCone,...
                'truncationTau',tau,...
                'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
                'truncationCircleRadius',effectiveRadii/tau);
        end

        function [optimalVelocity] = strategy_minimumDifference(desiredVelocity,escapeVelocities)
            inputDim = numel(desiredVelocity);
            searchMatrix = vertcat(escapeVelocities,zeros(1,size(escapeVelocities,2)));
            for i = 1:size(searchMatrix,2)
                searchMatrix(inputDim+1,i) = norm(desiredVelocity - searchMatrix(1:inputDim,i));
            end
            [~,minIndex] = min(searchMatrix(inputDim+1,:),[],2);
            optimalVelocity = searchMatrix(1:inputDim,minIndex);
            if isempty(optimalVelocity)
                warning('No viable velocities found in search matrix.');
                optimalVelocity = zeros(inputDim,1);
            end
        end

        function [cubePoints] = GetFeasabilityGrid(minVector,maxVector,pointDensity)
            dimensionality = numel(minVector);
            numPoints = pointDensity^dimensionality;
            dimMultipliers = zeros(dimensionality,1);
            dimensionalDistribution = zeros(dimensionality,pointDensity);
            for dim = 1:dimensionality
                dimensionalDistribution(dim,:) = linspace(minVector(dim),maxVector(dim),pointDensity);
                dimMultipliers(dim) = numPoints/pointDensity^dim;
            end
            dimMultipliers = fliplr(dimMultipliers);

            cubePoints = zeros(dimensionality,size(dimensionalDistribution,2)^dimensionality);
            for dim = 1:dimensionality
                insertVector = [];
                for i = 1:pointDensity
                    distributionIter = repmat(dimensionalDistribution(dim,i),1,dimMultipliers(dim));
                    insertVector = horzcat(insertVector,distributionIter);
                end
                no_copies = numPoints/size(insertVector,2);
                cubePoints(dim,:) = repmat(insertVector,1,no_copies);
            end
        end

        function [p_inter,isSuccessful] = twoRayIntersection2D(P1,dP1,P2,dP2)
            assert(numel(P1) == 2,'Input must be 2D');
            assert(numel(P2) == 2,'Input must be 2D');
            isSuccessful = false;
            p_inter = NaN(2,1);
            div = dP1(2)*dP2(1) - dP1(1)*dP2(2);
            if div == 0
                disp('Lines are parallel');
                return
            end
            mua = (dP2(1)*(P2(2) - P1(2)) + dP2(2)*(P1(1) - P2(1))) / div;
            mub = (dP1(1)*(P2(2) - P1(2)) + dP1(2)*(P1(1) - P2(1))) / div;
            if mua < 0 || mub < 0
                return
            end
            p_inter = P1 + mua*dP1;
            isSuccessful = true;
        end

        function [p_inter,isSuccessful] = findAny2DIntersection(P1,dP1,P2,dP2)
            assert(numel(P1) == 2,'Input must be 2D');
            assert(numel(P2) == 2,'Input must be 2D');
            isSuccessful = false;
            p_inter = NaN(2,1);
            div = dP1(2)*dP2(1) - dP1(1)*dP2(2);
            if div == 0
                return
            end
            mua = (dP2(1)*(P2(2) - P1(2)) + dP2(2)*(P1(1) - P2(1))) / div;
            p_inter = P1 + mua*dP1;
            isSuccessful = true;
        end

        function [projectedPoint,isOnTheRay] = pointProjectionToRay(p,p0,v0)
            projectedPoint = v0*v0'/(v0'*v0)*(p - p0) + p0;
            isOnTheRay = v0'*(projectedPoint - p0) > 0;
        end

        function [flag] = isInsideVO(point,VO)
            flag = 0;
            VOtolerance = 1E-8;
            candidateVector = point - VO.apex;
            VOprojection = norm(candidateVector)*cos(VO.openAngle/2);
            candProjection = VO.axisUnit'*candidateVector;
            if (candProjection - VOprojection) > VOtolerance
                flag = 1;
            end
        end

        function [pa,pb,isSuccessful] = findAny3DIntersection(P1,dP1,P3,dP3)
            EPS = 1E-9;
            isSuccessful = false;
            pa = NaN(3,1);
            pb = NaN(3,1);
            p13 = P1 - P3;
            d1343 = p13(1)*dP3(1) + p13(2)*dP3(2) + p13(3)*dP3(3);
            d4321 = dP3(1)*dP1(1) + dP3(2)*dP1(2) + dP3(3)*dP1(3);
            d1321 = p13(1)*dP1(1) + p13(2)*dP1(2) + p13(3)*dP1(3);
            d4343 = dP3(1)*dP3(1) + dP3(2)*dP3(2) + dP3(3)*dP3(3);
            d2121 = dP1(1)*dP1(1) + dP1(2)*dP1(2) + dP1(3)*dP1(3);
            denom = d2121*d4343 - d4321*d4321;
            if abs(denom) < EPS
                return
            end
            numer = d1343*d4321 - d1321*d4343;
            mua = numer/denom;
            mub = (d1343 + d4321*mua)/d4343;
            pa = P1 + mua*dP1;
            pb = P3 + mub*dP3;
            isSuccessful = true;
        end

        function [Cone] = vectorCone(pointA,pointB,radialPoint,nodes,coneColour)
            if nargin < 4 || isempty(nodes)
                nodes = 10;
            end
            coneEdgeColour = 'r';
            coneAlpha = 0.1;

            axisVector = pointB - pointA;
            mod_axisVector = sqrt(sum(axisVector.^2));
            tangent = radialPoint - pointA;
            mod_tangent = sqrt(sum(tangent.^2));
            trueAB = (dot(tangent,axisVector)/mod_axisVector^2)*axisVector;
            mod_trueAB = sqrt(sum(trueAB.^2));
            mod_radius = sqrt(mod_tangent^2 - mod_trueAB^2);

            t = linspace(0,2*pi,nodes)';
            xa2 = zeros(length(t),1);
            xa3 = zeros(size(xa2));
            xb2 = mod_radius*cos(t);
            xb3 = mod_radius*sin(t);
            x1 = [0 mod_trueAB];
            xx1 = repmat(x1,length(xa2),1);
            xx2 = [xa2 xb2];
            xx3 = [xa3 xb3];

            Cone = mesh(gca,real(xx1),real(xx2),real(xx3));
            unit_Vx = [1 0 0];
            angle_X1X2 = acos(dot(unit_Vx,axisVector)/(norm(unit_Vx)*mod_axisVector))*180/pi;
            axis_rot = cross([1 0 0],axisVector);
            if angle_X1X2 ~= 0
                rotate(Cone,axis_rot,angle_X1X2,[0 0 0])
            end
            set(Cone,'XData',get(Cone,'XData') + pointA(1))
            set(Cone,'YData',get(Cone,'YData') + pointA(2))
            set(Cone,'ZData',get(Cone,'ZData') + pointA(3))
            set(Cone,'AmbientStrength',1,...
                'FaceColor',coneColour,...
                'FaceLighting','gouraud',...
                'FaceAlpha',coneAlpha,...
                'EdgeColor',coneEdgeColour,...
                'EdgeAlpha',0);
        end
    end
end
