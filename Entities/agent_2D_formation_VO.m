
% Author: James A. Douthwaite

classdef agent_2D_formation_VO < agent_2D_VO & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTITIES UNIQUE TO AGENT CLASS
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_2D_formation_VO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            this = this@agent_2D_VO(varargin);
            
            this.maxSamples = 1;
            
            % //////////////////// SENSOR PARAMETERS //////////////////////
%             [obj.SENSORS] = obj.GetDefaultSensorParameters();       % Default sensing
            [this.SENSORS] = this.GetCustomSensorParameters();       % Experimental sensing
            % /////////////////////////////////////////////////////////////

            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % obj      - The agent object
            % TIME     - The current time structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated object
                        
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if this.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + this.objectID);
                hold on; grid on;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet,~] = this.GetAgentUpdate(ENV,varargin{1}); 

            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only.
            L = 0;
            if ~isempty(agentSet) 
                % PASS AGENT SET TO FORMATION CONTROLLER
                [heading,speed,L] = this.formationControl_distance(agentSet);          % Get the force vector
            
                desiredVelocity = heading*speed;
            end
           
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
%             avoidanceSet = [obstacleSet;agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  
%             if ~isempty(avoidanceSet) && avoidanceEnabled
%                 algorithm_indicator = 1;
%                 % GET THE UPDATED DESIRED VELOCITY
%                 [desiredHeadingVector,desiredSpeed] = this.GetAvoidanceCorrection(...
%                     dt,...
%                     desiredVelocity,...
%                     visualiseProblem);
%                 desiredVelocity = desiredHeadingVector*desiredSpeed;
%             end
            algorithm_dt = toc(algorithm_start); 
            
            % //////////////////// AGENT CONTROL //////////////////////////     
            [this] = this.Controller(ENV.dt,desiredVelocity);
                       
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            this.DATA.inputNames = {'$v_x$ (m/s)','$v_y$ (m/s)','$\dot{\psi}$ (rad/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(4:6);         % Record the control inputs
            this.DATA.lypanov(ENV.currentStep) = L;                        % Record the lypanov value            
        end
    end
end