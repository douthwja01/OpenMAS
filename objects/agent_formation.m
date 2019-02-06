%% FORMATION CONTROL AGENT TOOL SET (agent_formation.m) %%%%%%%%%%%%%%%%%%%
% This agent is designed to implement the formation control laws used in
% our current coordination studies. The agent has the ability to
% design and inact a control vector that drives the agents towards a
% formation defined by an adjacency matrix.

% Author: James A. Douthwaite & Dr Shiyu Zhao

classdef agent_formation < agent
    %%  INITIALISE THE PROPERTIES OF THE FORMATION CONTROLLER
    properties
        % FORMATION PARAMETERS
        adjacencyMatrix; % Example adjacency matrix
    end
    %%  CLASS METHODS
    methods
        % CONSTRUCTION METHOD
        function obj = agent_formation(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent(varargin);                                           % Get the supercalss
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
    end
    % ///////////// SHIYU'S FORMATION CONTROL TECHNIQUES //////////////////
    methods (Access = public)
        % BEARING BASED FORMATION CONTROL [INCOMPLETE]
        function [vi,V] = formationControl_bearing(obj,observedObjects)
            % This function calculates the formation control vector to
            % bring about the desired bearings between the given objects.
            
            % Author: Dr Shiyu Zhao
            
            % CONSTANT PARAMETERS
            objectNumber = numel(observedObjects);
            
            %             pi = obj.localState(1:numel(observedObjects(1).position));
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            for j = 1:objectNumber
                % SECOND AGENT INDEX
                objectID_B = observedObjects(j).objectID;
                pij = observedObjects(j).position;
                
                
                % ACCUMULATE THE BEARING FEEDBACK
                vi = vi + Pij*(pij);
                % GET THE LYAPANOV VALUE
                V = V + (norm(Pij*(pij)))^2;
            end
            % CONDITION THE OUTPUT VECTOR
            [vi] = obj.conditionControlVector(vi);
        end
        % DISTANCE BASE FORMATION CONTROL [COMPLETE]
        function [vi,V] = formationControl_distance(obj,observedObjects)
            % This function calculates the formation control vector to
            % to bring about the desired object separation, given the
            % agents current knowledge of the surrounding agents.
            
            % Author: Dr Shiyu Zhao
            
            % Under this implementation, the adjacency matrix is indexed by
            % agent objectID's. This means agent 2 will only calculate the
            % contributions from the observed agents 1 and 2. The resulting
            % adjacency matrix is of the form:
            % adj = [0   d2  d3;
            %        d2   0  d1;
            %        d3  d1   0];
            
            % CONFIRM ADJACENCY MATRIX
            if ~isprop(obj,'adjacencyMatrix') || isempty(obj.adjacencyMatrix)
                error('Agent is missing (or has not been assigned) an adjacency matrix');
            end
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedObjects);
            pi = obj.localState(1:numel(observedObjects(1).position));
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            for j = 1:agentNumber
                % SECOND AGENT INDEX
                objectID_B = observedObjects(j).objectID;
                pj = pi + observedObjects(j).position;
                % separation BASED CONFIRMATION
                ell_ij = obj.adjacencyMatrix(obj.objectID,objectID_B);
                % SUM THE CONTRIBUTION
                vi = vi + (norm(pi - pj)^2 - ell_ij^2)*(pj - pi);
                % GET THE LYAPANOV VALUE
                V = V + (norm(pi-pj)^2-ell_ij^2)^2;
            end
            % CONDITION THE OUTPUT VECTOR
            [vi] = obj.conditionControlVector(vi);
        end
        % RELATIVE POSITION BASED FORMATION CONTROLLER [INCOMPLETE]
        function [vi,V] = formationControl_relativePosition(obj,observedObjects)
            % This function calculates the formation control vector to
            % to bring about the desired object position, given the
            % agents current knowledge of the surrounding agents.
            
            % Author: Dr Shiyu Zhao
            
            % Under this implementation, the adjacency matrix is indexed by
            % agent objectID's. This means agent 2 will only calculate the
            % contributions from the observed agents 1 and 2.
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedObjects);
            
            %             pi = obj.localState(1:numel(observedObjects(1).position));
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            for j = 1:agentNumber
                % SECOND AGENT INDEX
                objectID_B = observedObjects(j).objectID;
                % THE AJACENCY MATRIX COMPONENT
                scale_ij = obj.adjacencyMatrix(obj.objectID,objectID_B);
                
                % RELATIVE POSITION
                pij = observedObjects(j).position;                          % The observed position is already relative
                
                % THE DESIRED RELATIVE POSITION OF AGENT j TO i
                pij_star = obj.relativePositionMatrix(obj.objectID,objectID_B);
                
                % ACCUMULATIVE FEEDBACK FROM THE SET OF RELATIVE POSITIONS
                vi = vi + scale_ij*(pij - pij_star);
                
                % GET THE LYAPANOV VALUE
                V = V + norm(scale_ij*(pij - pij_star))^2;
            end
            % CONDITION THE OUTPUT VECTOR
            [vi] = obj.conditionControlVector(vi);
        end
        % CONDITION OUTPUT VECTOR
        function [vi] = conditionControlVector(obj,vi)
            % Ensures the command vector is achievable and lies within the
            % constraints specified by 'obj.maxSpeed'. Works for both 2D
            % and 3D velocities.
            
            % DESIGN SUITABLE OUTPUT VECTOR
            norm_vi = norm(vi);
            if norm_vi == 0
                unit_vi = [1;zeros(numel(vi)-1,1)];
            else
                unit_vi = vi/norm_vi;
            end
            
            if isprop(obj,'nominalSpeed')
                norm_vi = obj.boundValue(norm_vi,-obj.nominalSpeed,obj.nominalSpeed);
                vi = unit_vi*norm_vi;
            else
                error('Nominal speed not defined');
            end
        end
    end
    % ///////////////// PRIMATIVE FLOCKING ALGORITHIMS
    methods
        % THE FORMATION CONTROL PRINCIPLES FOR THE 'BOIDS' MODEL
        function [v_boids] = formationControl_boids(obj,targetWaypoint,neighbours,weights)
            % This function computes the formation control heading based
            % upon the boids four principle rules.
            
            if nargin < 4
                weights = ones(4,1);    % Default to balances contributions
            end
            % APPLY THE RULE SET
            [v_sep] = obj.separationRule(neighbours);
            [v_ali] = obj.alignmentRule(neighbours);
            [v_coh] = obj.cohesionRule(neighbours);
            [v_mig] = obj.migrationRule(targetWaypoint);
            % CALCULATE THE WEIGHTED CONTRIBUTIONS
            v_boids = weights(1)*v_sep + weights(2)*v_ali + weights(3)*v_coh + weights(4)*v_mig;
            % RENORMALISE
            v_boids = OMAS_geometry.unit(v_boids);
            
%             v_boids = obj.VIRTUAL.R*v_boids
            
        end
    end
    
    methods (Static)
        % THE SEPARATION RULE
        function [v_sep] = separationRule(neighbours)
            % We assume the position of the neighbours are measured
            % relatively. We construct a vector that is the sum of
            %             v_sep = zeros(3,1);
            v_sep = [1;0;0];
            for n = 1:numel(neighbours)
                pn = neighbours(n).position;
                % Get the (-ve) separation vector
                if sum(abs(pn)) == 0
                    pn = [1;0;0]*1E-5;
                end
                vp = -pn;
                % Normalise
                vn = OMAS_geometry.unit(vp);
                % Scale by the proximity
                vn = vn*norm(vn);
                % Combine for net influence vector
                v_sep = v_sep + vn;
            end
            
        end
        % THE ALIGNMENT RULE
        function [v_ali] = alignmentRule(neighbours)
            % This rule urges the robots to move in the same direction as
            % its neighbours. (i.e heading matching)
            %             v_ali = zeros(3,1);
            v_ali = [1;0;0];
            for n = 1:numel(neighbours)
                % THE RELATIVE VELOCITY VECTORS
                v_ali = v_ali + neighbours(n).velocity;
            end
            v_ali = v_ali/numel(v_ali);
        end
        % THE COHESION RULE
        function [v_coh] = cohesionRule(neighbours)
            % This rule will try to move the agent towards the center of
            % mass of the other robots to form a group/swarm.
            %             v_coh = zeros(3,1);
            v_coh = [1;0;0];
            for n = 1:numel(neighbours)
                % SUM THE VECTOR POSITIONS
                v_coh = v_coh + neighbours(n).position;
            end
            % THE MEAN POSITION
            v_coh = v_coh/numel(neighbours);
        end
        % THE MIGRATION RULE
        function [v_mig] = migrationRule(targetWaypoint)
            % This rule aims to move the swarm/individual towards a
            % designated goal location.
            if ~isstruct(targetWaypoint)
                v_mig = [1;0;0];
                return
            end
            % THE RELATIVE POSITION VECTOR
            v_mig = targetWaypoint.position;
        end
    end
% methods
%     % THE SEPARATION RULE
%     function [v_sep] = separationRule(obj,neighbours)
%         % We assume the position of the neighbours are measured
%         % relatively. We construct a vector that is the sum of
%         %             v_sep = zeros(3,1);
%         v_sep = [1;0;0];
%         for n = 1:numel(neighbours)
%             pn = obj.VIRTUAL.globalPosition - neighbours(n).globalPosition;
%             % Get the (-ve) separation vector
%             if sum(abs(pn)) == 0
%                 pn = [1;0;0]*1E-5;
%             end
%             % Normalise
%             vn = OMAS_geometry.unit(pn);
%             % Scale by the proximity
%             vn = vn*norm(pn);
%             % Combine for net influence vector
%             v_sep = v_sep - vn;
%         end
%     end
%     % THE ALIGNMENT RULE
%     function [v_ali] = alignmentRule(obj,neighbours)
%         % This rule urges the robots to move in the same direction as
%         % its neighbours. (i.e heading matching)
%         %             v_ali = zeros(3,1);
%         v_ali = [1;0;0];
%         for n = 1:numel(neighbours)
%             % THE RELATIVE VELOCITY VECTORS
%             v_ali = v_ali + neighbours(n).globalVelocity;
%         end
%         % The mean alignment
%         v_ali = v_ali/numel(neighbours);
%     end
%     % THE COHESION RULE
%     function [v_coh] = cohesionRule(obj,neighbours)
%         % This rule will try to move the agent towards the center of
%         % mass of the other robots to form a group/swarm.
%         %             v_coh = zeros(3,1);
%         v_coh = [1;0;0];
%         for n = 1:numel(neighbours)
%             % SUM THE VECTOR POSITIONS
%             v_coh = v_coh + neighbours(n).globalPosition;
%         end
%         % THE MEAN POSITION
%         v_coh = v_coh/numel(neighbours);
%     end
%     % THE MIGRATION RULE
%     function [v_mig] = migrationRule(obj,targetWaypoint)
%         % This rule aims to move the swarm/individual towards a
%         % designated goal location.
%         if ~isstruct(targetWaypoint)
%             v_mig = [1;0;0];
%             return
%         end
%         % THE RELATIVE POSITION VECTOR
%         v_mig = targetWaypoint.globalPosition - obj.VIRTUAL.globalPosition;
%     end
% end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]