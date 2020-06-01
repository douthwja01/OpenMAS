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
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]

