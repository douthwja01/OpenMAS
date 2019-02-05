
% Author: James A. Douthwaite

classdef agent_2D_formation_VO < agent_2D_VO & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTITIES UNIQUE TO AGENT CLASS
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_2D_formation_VO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_2D_VO(varargin);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);   
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
        % This is necessary to force the agent_2D_VO/RVO/HRVO/RVO2 to
        % consider the formation control desired velocity vector.
        function [obj] = main(obj,TIME,varargin)
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
            dt = TIME.dt;
            desiredSpeed = obj.nominalSpeed;
            desiredHeadingVector = [1;0;0];
            desiredVelocity = desiredHeadingVector*desiredSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.getAgentUpdate(varargin{1});       % IDEAL INFORMATION UPDATE

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
                [desiredHeadingVector,desiredSpeed] = obj.getAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start); 
            
            % //////////////////// AGENT CONTROL //////////////////////////            
            % APPLY SPEED CONSTRAINT
            desiredSpeed = norm(desiredVelocity);
            if desiredSpeed > obj.maxSpeed
                desiredHeadingVector = desiredVelocity/norm(desiredVelocity);
                desiredVelocity = desiredHeadingVector*obj.maxSpeed;
            elseif desiredSpeed == 0
                desiredHeadingVector = [1;0];
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dHeading] = obj.getVectorHeadingAngles([1;0;0],[desiredHeadingVector;0]); % Relative heading angles   
            omega = -dHeading/dt;
            desiredVelocity = [desiredSpeed;0];
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:3),desiredVelocity,omega);
            obj.localState(1:3) = obj.localState(1:3) + dt*dX;
            obj.localState(4:6) = dX;
            
            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES /////////////// 
            [obj] = obj.updateGlobalProperties_2DVelocities(dt,obj.localState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(TIME,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'dx (m/s)','dy (m/s)','Yaw Rate (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = obj.localState(4:6);         % Record the control inputs
            obj.DATA.lypanov(TIME.currentStep) = L;                               % Record the lypanov value
        end
    end
end