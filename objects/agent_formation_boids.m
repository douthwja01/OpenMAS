
% Author: James A. Douthwaite

classdef agent_formation_boids < agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        weights = [0.1;5;1;0.2];        % Weights for the boids behaviour
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_formation_boids(varargin)
            % Constructor
            this = this@agent_formation(varargin);
            
            % For clarity ->DYNAMICS
            this.v_max = 5;
            this.w_max = 2;
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % ///////////////////////////////////////////////////////////// 
        end
        % Setup
        % - The same as any 2D/3D object
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % this      - The agent object
            % TIME     - The current time structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % this      - The updated object
                        
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + this.objectID);
                hold on; grid on;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % DEFAULT BEHAVIOUR 
            desiredSpeed    = this.v_nominal;
            desiredHeadingVector = [1;0;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;

            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,~,agentSet] = this.GetAgentUpdate(ENV,varargin{1});         % IDEAL INFORMATION UPDATE

            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only.
            algorithm_start = tic; algorithm_indicator = 0; 
            if ~isempty(agentSet) 
                algorithm_indicator = 1;
                % Formation control routine
%                 [v_boids] = this.formationControl_boids(...
%                     this.targetWaypoint,...
%                     agentSet,...
%                     this.weights);

                [v_boids] = this.formationControl_boids(...
                    this.targetWaypoint,...
                    agentSet,...
                    this.weights);
                
                desiredVelocity = desiredSpeed*v_boids; % Define the velocity request
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm      
            
            % Pass the velocity to the controller
            this = this.Controller(ENV.dt,desiredVelocity);
                        
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            this.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = [this.localState(7);this.localState(4:6)];         % Record the control inputs
        end
        
    end
    
    %% BOIDS FLOCKING ALGORITHIMS
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