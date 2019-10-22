
% Author: James A. Douthwaite

classdef agent_2D_formation_VO < agent_2D_VO & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTITIES UNIQUE TO AGENT CLASS
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function obj = agent_2D_formation_VO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_2D_VO(varargin);
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % obj      - The agent object
            % TIME     - The current time structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated object
                        
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + obj.objectID);
                hold on; grid on;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % DEFAULT BEHAVIOUR 
            dt = ENV.dt;
            
            desiredSpeed = obj.nominalSpeed;
            desiredHeadingVector = [1;0;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE

            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only.
            L = 0;
            if ~isempty(agentSet) 
                % PASS AGENT SET TO FORMATION CONTROLLER
                [vi,L] = obj.formationControl_distance(agentSet);          % Get the force vector
                % DEFINE VELOCITY REQUEST
                desiredVelocity = vi;
            end
                         
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet;agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = obj.GetAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start); 
            
            % //////////////////// AGENT CONTROL //////////////////////////     
            [obj] = obj.controller(dt,desiredVelocity);
                       
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'dx (m/s)','dy (m/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = obj.localState(4:6);         % Record the control inputs
            obj.DATA.lypanov(ENV.currentStep) = L;                               % Record the lypanov value
        end
    end
end