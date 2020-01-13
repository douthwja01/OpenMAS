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
        function this = agent_formation(varargin)
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            this@agent(varargin);                                           % Get the supercalss
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
    end
    % ///////////// SHIYU'S FORMATION CONTROL TECHNIQUES //////////////////
    methods (Access = public)
        % BEARING BASED FORMATION CONTROL [INCOMPLETE]
        function [vi,V] = formationControl_bearing(this,observedObjects)
            % This function calculates the formation control vector to
            % bring about the desired bearings between the given objects.
            
            % Author: Dr Shiyu Zhao
            
            % CONSTANT PARAMETERS
            objectNumber = numel(observedObjects);
            
            % pi = this.localState(1:numel(observedObjects(1).position));
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
            [vi] = this.conditionControlVector(vi);
        end
        % DISTANCE BASE FORMATION CONTROL [COMPLETE]
        function [heading,speed,V] = formationControl_distance(this,observedObjects)
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
            if ~isprop(this,'adjacencyMatrix') || isempty(this.adjacencyMatrix)
                error('Agent is missing (or has not been assigned) an adjacency matrix');
            end
            
            if this.Is3D()
                pi = this.localState(1:3);
                vi = zeros(3,1);
            else
                pi = this.localState(1:2);
                vi = zeros(2,1);
            end
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedObjects);
            V = 0;
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            for j = 1:agentNumber
                % GET THE SECOND OBJECT DATA
                objectID_j = this.GetLastMeasurementFromStruct(observedObjects(j),'objectID');
                p_j = this.GetLastMeasurementFromStruct(observedObjects(j),'position');
                % SECOND AGENT INDEX
                pj = pi + p_j;
                % separation BASED CONFIRMATION
                ell_ij = this.adjacencyMatrix(this.objectID,objectID_j);
                % SUM THE CONTRIBUTION
                vi = vi + (norm(pi - pj)^2 - ell_ij^2)*(pj - pi);
                % GET THE LYAPANOV VALUE
                V = V + (norm(pi-pj)^2-ell_ij^2)^2;
            end
            % CONDITION THE OUTPUT VECTOR
            %[vi] = this.conditionControlVector(vi);
            speed   = norm(vi);
            heading = vi/speed;
        end
        % RELATIVE POSITION BASED FORMATION CONTROLLER [INCOMPLETE]
        function [vi,V] = formationControl_relativePosition(this,observedObjects)
            % This function calculates the formation control vector to
            % to bring about the desired object position, given the
            % agents current knowledge of the surrounding agents.
            
            % Author: Dr Shiyu Zhao
            
            % Under this implementation, the adjacency matrix is indexed by
            % agent objectID's. This means agent 2 will only calculate the
            % contributions from the observed agents 1 and 2.
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedObjects);
            
            %             pi = this.localState(1:numel(observedObjects(1).position));
            vi = zeros(size(observedObjects(1).position));
            V = 0;
            for j = 1:agentNumber
                % SECOND AGENT INDEX
                objectID_B = observedObjects(j).objectID;
                % THE AJACENCY MATRIX COMPONENT
                scale_ij = this.adjacencyMatrix(this.objectID,objectID_B);
                
                % RELATIVE POSITION
                pij = observedObjects(j).position;                          % The observed position is already relative
                
                % THE DESIRED RELATIVE POSITION OF AGENT j TO i
                pij_star = this.relativePositionMatrix(this.objectID,objectID_B);
                
                % ACCUMULATIVE FEEDBACK FROM THE SET OF RELATIVE POSITIONS
                vi = vi + scale_ij*(pij - pij_star);
                
                % GET THE LYAPANOV VALUE
                V = V + norm(scale_ij*(pij - pij_star))^2;
            end
            % CONDITION THE OUTPUT VECTOR
            [vi] = this.conditionControlVector(vi);
        end
        % CONDITION OUTPUT VECTOR
        function [vi] = conditionControlVector(this,vi)
            % Ensures the command vector is achievable and lies within the
            % constraints specified by 'this.maxSpeed'. Works for both 2D
            % and 3D velocities.
            
            % DESIGN SUITABLE OUTPUT VECTOR
            norm_vi = norm(vi);
            if norm_vi == 0
                unit_vi = [1;zeros(numel(vi)-1,1)];
            else
                unit_vi = vi/norm_vi;
            end
            % Bound the value between min and maximum velocities
            norm_vi = boundValue(norm_vi,-this.v_nominal,this.v_nominal);
            % Return the normalised velocity
            vi = unit_vi*norm_vi;
        end
    end
    % ///////////////// PRIMATIVE FLOCKING ALGORITHIMS
    methods
        % THE FORMATION CONTROL PRINCIPLES FOR THE 'BOIDS' MODEL
        function [v_boids] = formationControl_boids(this,targetWaypoint,neighbours,weights)
            % This function computes the formation control heading based
            % upon the boids four principle rules.
            
            if nargin < 4
                weights = ones(4,1);    % Default to balances contributions
            end
            
            p_j = []; v_j = [];
            for j = 1:numel(neighbours)
                p_j(j,:) = this.GetLastMeasurementFromStruct(neighbours(j),'position')';
                v_j(j,:) = this.GetLastMeasurementFromStruct(neighbours(j),'velocity')';
            end
            % The way-point position
            if ~isempty(targetWaypoint)
                p_wp = this.GetLastMeasurementFromStruct(targetWaypoint,'position');
            else
                p_wp = this.localState(1:3);
            end
            
            % APPLY THE RULE SET
            [v_sep] = this.separationRule(p_j);
            [v_ali] = this.alignmentRule(v_j);
            [v_coh] = this.cohesionRule(p_j);
            [v_mig] = this.migrationRule(p_wp);
            % Calculated the weighted vectors
            v_boids = weights(1)*v_sep + weights(2)*v_ali + weights(3)*v_coh + weights(4)*v_mig;
            % Renormalise the vector
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
            assert(IsColumn(position_wp),'Expecting a vector of positions [n x dim].');
            
            % THE RELATIVE POSITION VECTOR
            v_mig = position_wp;
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]

