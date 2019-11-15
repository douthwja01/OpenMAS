%% VELOCITY OBSTACLE/FORMATION CONTROL AGENT (agent_formation_VO.m) %%%%%%%
% This agent is a hybrid of the formation control element agent class and
% the Velocity Obstacle (VO) class. 

% Author: James A. Douthwaite

classdef agent_formation_VO < agent_VO & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        
    end
    %% ////////////////////// MAIN CLASS METHODS //////////////////////////
    methods 
        % Constructor
        function [this] = agent_formation_VO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            this = this@agent_VO(varargin);
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % ///////////////////////////////////////////////////////////// 
        end
        % Setup X = [x;x_dot]' 
        function [this] = setup(this,localXYZVelocity,localXYZrotations)
            % The state initialiser must be called 'initialise_localState'
            % and instead calls the 'initialise_3DVelocities' function in
            % this case. 
            
            % BUILD THE STATE VECTOR FOR A 3D SYSTEM WITH CONCATINATED VELOCITIES
            [this] = this.setup_3DVelocities(localXYZVelocity,localXYZrotations);
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
            
            % DEFAULT BEHAVIOUR 
            desiredHeadingVector = [1;0;0];
            desiredVelocity = desiredHeadingVector*this.v_nominal;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE

            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only. 
            L = 0;
            if ~isempty(agentSet) 
                % PASS AGENT SET TO FORMATION CONTROLLER
                [desiredVelocity,L] = this.formationControl_distance(agentSet);      % Get the force vector
            end
            
            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            % Modify the desired velocity with the augmented avoidance velocity.
            avoidanceSet = [obstacleSet;agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredHeadingVector,desiredSpeed] = this.GetAvoidanceCorrection(desiredVelocity,avoidanceSet,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm 
                        
            % /////// COMPUTE STATE CHANGE FROM CONTROL INPUTS ////////////
            [this] = this.Controller(ENV.dt,desiredVelocity);
                        
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this = this.writeAgentData(ENV,algorithm_indicator,algorithm_dt);
            this.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = [this.localState(7);this.localState(4:6)];         % Record the control inputs
            this.DATA.lypanov(ENV.currentStep) = L;                               % Record the lypanov value
        end
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]