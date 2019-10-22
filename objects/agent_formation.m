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
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % Constructor
        function obj = agent_formation(varargin)
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent(varargin);                                           % Get the supercalss
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
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
            
            % pi = obj.localState(1:numel(observedObjects(1).position));
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
            
            if obj.VIRTUAL.is3D
                pi = obj.localState(1:3);
                vi = zeros(3,1);
            else
                pi = obj.localState(1:2);
                vi = zeros(2,1);
            end
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedObjects);
            V = 0;
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            for j = 1:agentNumber
                % GET THE SECOND OBJECT DATA
                objectID_j = obj.GetLastMeasurementFromStruct(observedObjects(j),'objectID');
                p_j = obj.GetLastMeasurementFromStruct(observedObjects(j),'positions');
                % SECOND AGENT INDEX
                pj = pi + p_j;
                % separation BASED CONFIRMATION
                ell_ij = obj.adjacencyMatrix(obj.objectID,objectID_j);
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
            p_j = []; v_j = [];
            for j = 1:numel(neighbours)
                p_j(j,:) = obj.GetLastMeasurementFromStruct(neighbours(j),'position')';
                v_j(j,:) = obj.GetLastMeasurementFromStruct(neighbours(j),'velocity')';
            end
            % The way-point position
            p_wp = obj.GetLastMeasurementFromStruct(targetWaypoint,'position');
            
            % APPLY THE RULE SET
            [v_sep] = obj.separationRule(p_j);
            [v_ali] = obj.alignmentRule(v_j);
            [v_coh] = obj.cohesionRule(p_j);
            [v_mig] = obj.migrationRule(p_wp);
            
            % CALCULATE THE WEIGHTED CONTRIBUTIONS
            v_boids = weights(1)*v_sep + weights(2)*v_ali + weights(3)*v_coh + weights(4)*v_mig;
            
            % RENORMALISE
            v_boids = v_boids/norm(v_boids);
        end
    end
    
    methods (Static)
        % THE SEPARATION RULE
        function [v_sep] = separationRule(positions)
            % We assume the position of the neighbours are measured
            % relatively. We construct a vector that is the sum of

            % Input sanity check
            assert(isnumeric(positions),'Expecting a vector of positions [n x dim].');
            
            v_sep = [1;0;0];
            for n = 1:size(positions,1)
                pn = positions(n,:)';
                % Get the (-ve) separation vector
                if sum(abs(pn)) == 0
                    pn = [1;0;0]*1E-5;
                end
                vp = -pn;
                % Normalise
                vn = unit(vp);
                % Scale by the proximity
                vn = vn*norm(vn);
                % Combine for net influence vector
                v_sep = v_sep + vn;
            end
        end
        % THE ALIGNMENT RULE
        function [v_ali] = alignmentRule(velocities)
            % This rule urges the robots to move in the same direction as
            % its neighbours. (i.e heading matching)
            
            % Input sanity check
            assert(isnumeric(velocities),'Expecting a vector of positions [n x dim].');
            
            v_ali = [1;0;0];
            for n = 1:size(velocities,1)
                % THE RELATIVE VELOCITY VECTORS
                v_ali = v_ali + velocities(n,:)';
            end
            v_ali = v_ali/numel(v_ali);
        end
        % THE COHESION RULE
        function [v_coh] = cohesionRule(positions)
            % This rule will try to move the agent towards the center of
            % mass of the other robots to form a group/swarm.
            
            % Input sanity check
            assert(isnumeric(positions),'Expecting a vector of positions [n x dim].');
            
            v_coh = [1;0;0];
            for n = 1:size(positions,1)
                % SUM THE VECTOR POSITIONS
                v_coh = v_coh + positions(n,:)';
            end
            % THE MEAN POSITION
            v_coh = v_coh/size(positions,1);
        end
        % THE MIGRATION RULE
        function [v_mig] = migrationRule(position_wp)
            % This rule aims to move the swarm/individual towards a
            % designated goal location.
            
            % Input sanity check
            assert(isColumn(position_wp),'Expecting a vector of positions [n x dim].');
            
            % THE RELATIVE POSITION VECTOR
            v_mig = position_wp;
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]

