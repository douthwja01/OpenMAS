
% Author: James A. Douthwaite

classdef agent_formation_boids < agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        weights = [0.1;0.1;1;1];        % Weights for the boids behaviour
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_formation_boids(varargin)
            % Constructor
            this = this@agent_formation(varargin);
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
                [v_boids] = this.formationControl_boids(this.targetWaypoint,agentSet,this.weights);
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
end