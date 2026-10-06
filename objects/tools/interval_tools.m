%% INTERVAL TOOLS (interval_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared interval arithmetic and sensing used by interval agents.

classdef interval_tools
    methods (Static)
        function [SENSORS] = GetCustomSensorParameters()
            SENSORS = struct();
            SENSORS.range = 10;
            SENSORS.sigma_position = 0.5;
            SENSORS.sigma_velocity = 0.1;
            SENSORS.sigma_rangeFinder = 0.1;
            SENSORS.sigma_camera = 5.208E-5;
            SENSORS.sampleFrequency = inf;
            SENSORS.confidenceAssumption = 3;
        end

        function [SENSORS] = GetDefaultSensorParameters()
            SENSORS = struct();
            SENSORS.range = inf;
            SENSORS.sigma_position = 0.0;
            SENSORS.sigma_velocity = 0.0;
            SENSORS.sigma_rangeFinder = 0.0;
            SENSORS.sigma_camera = 0;
            SENSORS.sampleFrequency = inf;
            SENSORS.confidenceAssumption = 0;
        end

        function [observedObject] = SensorModel(agent,dt,observedObject)
            if isempty(observedObject)
                return
            end

            [psi_j,theta_j,alpha_j] = agent.GetCameraMeasurements(observedObject);
            azimuthBox   = midrad(psi_j,  agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_camera);
            elevationBox = midrad(theta_j,agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_camera);
            alphaBox     = midrad(alpha_j,agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_camera);

            d_j = agent.GetRangeFinderMeasurements(observedObject);
            rangeBox = midrad(d_j,agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_rangeFinder);
            rangeBox = intersect(infsup(0,inf),rangeBox);
            radiusBox = agent.GetRadiusFromAngularWidth(rangeBox,alphaBox);

            observedObject.range     = rangeBox;
            observedObject.heading   = azimuthBox;
            observedObject.elevation = elevationBox;
            observedObject.width     = alphaBox;
            observedObject.radius    = radiusBox;

            positionBox = GetCartesianFromSpherical(rangeBox,azimuthBox,elevationBox);
            if ~agent.Is3D
                positionBox = positionBox(1:2,1);
            end
            observedObject.position = positionBox;
            observedObject.velocity = observedObject.velocity + midrad(0,3*agent.SENSORS.sigma_velocity);
        end

        function [p_i,v_i,r_i] = GetAgentMeasurements(agent)
            assert(isstruct(agent.SENSORS) && ~isempty(agent.SENSORS),'The SENSOR structure is absent.');
            assert(isnumeric(agent.localState),'The state of the object is invalid.');
            assert(~any(isnan(agent.localState)),'The state contains NaNs');

            if agent.Is3D
                positionIndices = 1:3;
                velocityIndices = 7:9;
            else
                positionIndices = 1:2;
                velocityIndices = 4:5;
            end

            p_uncertainty = agent.SENSORS.sigma_position*randn(numel(positionIndices),1);
            v_uncertainty = agent.SENSORS.sigma_velocity*randn(numel(velocityIndices),1);
            p_i = agent.localState(positionIndices,1) + midrad(p_uncertainty,agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_position);
            v_i = agent.localState(velocityIndices,1) + midrad(v_uncertainty,agent.SENSORS.confidenceAssumption*agent.SENSORS.sigma_velocity);
            r_i = agent.radius;
        end

        function [headingVector] = GetTargetHeading(agent,targetObject)
            if nargin > 1
                targetPosition = agent.GetLastMeasurement(targetObject.objectID,'position');
            elseif ~isempty(agent.targetWaypoint)
                targetPosition = agent.targetWaypoint.position(:,agent.targetWaypoint.sampleNum);
            elseif agent.Is3D()
                targetPosition = [1;0;0];
            else
                targetPosition = [1;0];
            end
            headingVector = interval_tools.iheading(targetPosition);
        end

        function [targetLogical] = GetTargetCondition(agent)
            currentPosition = agent.targetWaypoint.position(:,agent.targetWaypoint.sampleNum);
            currentRadius   = agent.targetWaypoint.radius(:,agent.targetWaypoint.sampleNum);
            targetLogical = 0 > interval_tools.inorm(mid(currentPosition)) - (mid(currentRadius) + agent.radius);
        end

        function [agent] = UpdateMemoryOrderByField(agent,field)
            assert(ischar(field),'Memory sort method must be a string.');
            assert(isfield(agent.MEMORY,field),'Field must belong to the memory structure.');
            if isintval([agent.MEMORY.(field)])
                [~,ind] = sort(mid([agent.MEMORY.(field)]),2,'descend');
            else
                [~,ind] = sort([agent.MEMORY.(field)],2,'descend');
            end
            agent.MEMORY = agent.MEMORY(ind);
        end

        function [priority_j] = GetObjectPriority(agent,objectID)
            logicalIDIndex = [agent.MEMORY(:).objectID] == objectID;
            priority_j = 1/mid(agent.MEMORY(logicalIDIndex).range(agent.MEMORY(logicalIDIndex).sampleNum));
        end

        function [memStruct] = CreateMEMORY(agent,horizonSteps)
            IntLab();
            if nargin < 2 || isempty(horizonSteps)
                horizonSteps = 10;
            end
            if agent.Is3D
                dim = 3;
            else
                dim = 2;
            end

            memStruct = struct();
            memStruct.name = '';
            memStruct.objectID = uint8(0);
            memStruct.type = OMAS_objectType.misc;
            memStruct.sampleNum = uint8(1);
            memStruct.time      = circularBuffer(NaN(1,horizonSteps));
            memStruct.position  = circularBuffer(intval(NaN(dim,horizonSteps)));
            memStruct.velocity  = circularBuffer(intval(NaN(dim,horizonSteps)));
            memStruct.radius    = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.range     = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.heading   = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.elevation = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.width     = circularBuffer(intval(NaN(1,horizonSteps)));
            memStruct.geometry = struct('vertices',[],'faces',[],'normals',[],'centroid',[]);
            memStruct.priority = [];
        end

        function [p] = linKin_position(dt,p0,v0)
            p = p0 + v0*dt;
        end

        function [v] = linKin_velocity(dt,p0,p1)
            v = (p1-p0)/dt;
        end

        function [lambda,theta] = GetVectorHeadingAngles(V,U)
            V = V/norm(V);
            U = U/norm(U);
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];
            rotationAxis = cross(Vh,Uh);
            lambda = sign(mid(rotationAxis(3)))*acos(dot(Vh,Uh)/norm(Vh));
            if numel(V) == 3 && numel(U) == 3
                theta = atan2(U(3),norm(Uh));
            else
                theta = 0;
            end
        end

        function [v_unit] = iheading(v)
            if ~isintval(v)
                v_unit = v/norm(v);
                return
            end
            v_sup = sup(v);
            v_inf = inf(v);
            v_sup_unit = v_sup/norm(v_sup);
            v_inf_unit = v_inf/norm(v_inf);
            v_unit = infsup(-1,1);
            for i = 1:size(v,1)
                if v_inf_unit(i) > v_sup_unit(i)
                    v_unit(i,1) = infsup(v_sup_unit(i),v_inf_unit(i));
                else
                    v_unit(i,1) = infsup(v_inf_unit(i),v_sup_unit(i));
                end
            end
        end

        function [v_int] = iunit(v)
            vnorm = interval_tools.inorm(v);
            if inf(vnorm) == 0
                vnorm = infsup(1E-12,sup(vnorm));
            end
            v_int = interval_tools.idivide(v,vnorm);
        end

        function [v_int] = idot(v_a,v_b)
            v_int = 0;
            for i = 1:size(v_a,1)
                v_int = v_int + v_a(i)*v_b(i);
            end
        end

        function [v_det] = idet(v_a,v_b)
            v_det = v_a(1)*v_b(2) - v_b(1)*v_a(2);
        end

        function [v_int] = icross(v_a,v_b)
            v_int = [...
                v_a(2)*v_b(3) - v_a(3).*v_b(2);
                v_a(3)*v_b(1) - v_a(1).*v_b(3);
                v_a(1)*v_b(2) - v_a(2).*v_b(1)];
        end

        function [n_int] = inorm(v)
            if ~isintval(v)
                n_int = norm(v);
                return
            end
            if sum(iszero(v)) == length(v)
                n_int = intval(1E-8);
                return
            end

            infVal = inf(v).^2;
            supVal = sup(v).^2;
            sqrMatrix = zeros(size(infVal,1),2);
            for ind = 1:length(infVal)
                if infVal(ind) > supVal(ind)
                    sqrMatrix(ind,:) = [supVal(ind),infVal(ind)];
                else
                    sqrMatrix(ind,:) = [infVal(ind),supVal(ind)];
                end
            end
            sqrtMatrix = sqrt(sum(sqrMatrix,1));
            if sqrtMatrix(1) == 0
                n_int = infsup(1E-8,sqrtMatrix(2));
            else
                n_int = infsup(sqrtMatrix(1),sqrtMatrix(2));
            end
        end

        function [v_int] = idivide(v_a,v_b)
            v_int = v_a./v_b;
            ii = find(inf(v_b).*sup(v_b) <= 0);
            infinity = 999999999999;
            for i = ii
                if inf(v_b(i)) < 0 && sup(v_b(i)) > 0
                    v_int(i) = infsup(-infinity,infinity);
                elseif inf(v_b(i)) == 0 && sup(v_b(i)) ~= 0
                    v_int(i) = v_a(i)*infsup(1/sup(v_b(i)),infinity);
                elseif inf(v_b(i)) ~= 0 && sup(v_b(i)) == 0
                    v_int(i) = v_a(i)*infsup(-infinity,inf(v_b(i)));
                end
            end
        end
    end
end
