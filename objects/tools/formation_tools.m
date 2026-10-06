%% FORMATION TOOLS (formation_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared formation and boids laws used by 2D, 3D and vehicle agents.

classdef formation_tools
    methods (Static)
        function [agent] = step3D(agent,ENV,observations)
            visualiseProblem = 0;
            desiredVelocity = [1;0;0]*agent.v_nominal;

            [agent,obstacleSet,agentSet] = agent.GetAgentUpdate(ENV,observations);

            L = 0;
            if ~isempty(agentSet)
                [desiredVelocity,L] = formation_tools.distance(agent,agentSet);
            end

            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic;
            algorithm_indicator = 0;
            if ~isempty(avoidanceSet)
                algorithm_indicator = 1;
                [desiredHeadingVector,desiredSpeed] = agent.GetAvoidanceCorrection(desiredVelocity,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);

            agent = agent.Controller(ENV.dt,desiredVelocity);
            agent = agent.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            agent.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            agent.DATA.inputs(1:length(agent.DATA.inputNames),ENV.currentStep) = [agent.localState(7);agent.localState(4:6)];
            agent.DATA.lypanov(ENV.currentStep) = L;
        end

        function [agent] = step2D(agent,ENV,observations)
            [agent,~,agentSet] = agent.GetAgentUpdate(ENV,observations);

            L = 0;
            desiredVelocity = [1;0]*agent.v_nominal;
            if ~isempty(agentSet)
                [heading,speed,L] = formation_tools.distance(agent,agentSet);
                desiredVelocity = heading*speed;
            end

            algorithm_start = tic;
            algorithm_indicator = 0;
            algorithm_dt = toc(algorithm_start);

            agent = agent.Controller(ENV.dt,desiredVelocity);
            agent = agent.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            agent.DATA.inputNames = {'$v_x$ (m/s)','$v_y$ (m/s)','$\dot{\psi}$ (rad/s)'};
            agent.DATA.inputs(1:length(agent.DATA.inputNames),ENV.currentStep) = agent.localState(4:6);
            agent.DATA.lypanov(ENV.currentStep) = L;
        end

        function [vi,V] = bearing(agent,observedObjects)
            objectNumber = numel(observedObjects);
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            for j = 1:objectNumber
                pij = observedObjects(j).position;
                vi = vi + Pij*(pij);
                V = V + (norm(Pij*(pij)))^2;
            end
            vi = formation_tools.condition(agent,vi);
        end

        function [heading,speed,V] = distance(agent,observedObjects)
            if ~isprop(agent,'adjacencyMatrix') || isempty(agent.adjacencyMatrix)
                error('Agent is missing (or has not been assigned) an adjacency matrix');
            end

            if agent.Is3D()
                pi = agent.localState(1:3);
                vi = zeros(3,1);
            else
                pi = agent.localState(1:2);
                vi = zeros(2,1);
            end

            V = 0;
            for j = 1:numel(observedObjects)
                objectID_j = agent.GetLastMeasurementFromStruct(observedObjects(j),'objectID');
                p_j = agent.GetLastMeasurementFromStruct(observedObjects(j),'position');
                pj = pi + p_j;
                ell_ij = agent.adjacencyMatrix(agent.objectID,objectID_j);
                vi = vi + (norm(pi - pj)^2 - ell_ij^2)*(pj - pi);
                V = V + (norm(pi-pj)^2 - ell_ij^2)^2;
            end
            speed = norm(vi);
            heading = vi/speed;
        end

        function [vi,V] = relativePosition(agent,observedObjects)
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            for j = 1:numel(observedObjects)
                objectID_B = observedObjects(j).objectID;
                scale_ij = agent.adjacencyMatrix(agent.objectID,objectID_B);
                pij = observedObjects(j).position;
                pij_star = agent.relativePositionMatrix(agent.objectID,objectID_B);
                vi = vi + scale_ij*(pij - pij_star);
                V = V + norm(scale_ij*(pij - pij_star))^2;
            end
            vi = formation_tools.condition(agent,vi);
        end

        function [vi] = condition(agent,vi)
            norm_vi = norm(vi);
            if norm_vi == 0
                unit_vi = [1;zeros(numel(vi)-1,1)];
            else
                unit_vi = vi/norm_vi;
            end
            norm_vi = boundValue(norm_vi,-agent.v_nominal,agent.v_nominal);
            vi = unit_vi*norm_vi;
        end

        function [v_sep] = separationRule(positions)
            assert(isnumeric(positions),'Expecting a vector of positions [n x dim].');
            v_sep = [1;0;0];
            for n = 1:size(positions,1)
                pn = positions(n,:)';
                if sum(abs(pn)) == 0
                    pn = [1;0;0]*1E-5;
                end
                vn = unit(-pn);
                v_sep = v_sep + vn*norm(vn);
            end
        end

        function [v_ali] = alignmentRule(velocities)
            assert(isnumeric(velocities),'Expecting a vector of positions [n x dim].');
            v_ali = [1;0;0];
            for n = 1:size(velocities,1)
                v_ali = v_ali + velocities(n,:)';
            end
            v_ali = v_ali/numel(v_ali);
        end

        function [v_coh] = cohesionRule(positions)
            assert(isnumeric(positions),'Expecting a vector of positions [n x dim].');
            v_coh = [1;0;0];
            for n = 1:size(positions,1)
                v_coh = v_coh + positions(n,:)';
            end
            v_coh = v_coh/size(positions,1);
        end

        function [v_mig] = migrationRule(position_wp)
            assert(IsColumn(position_wp),'Expecting a vector of positions [n x dim].');
            v_mig = position_wp;
        end
    end
end
