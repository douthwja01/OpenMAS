
% Author: James A. Douthwaite

classdef agent_formation_boids < agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        % PROPERTITIES UNIQUE TO AGENT CLASS
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_formation_boids(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_formation(varargin);
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
            [obj,~,agentSet] = obj.GetAgentUpdate(dt,varargin{1});         % IDEAL INFORMATION UPDATE

            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only.
            algorithm_start = tic; algorithm_indicator = 0; 
            if ~isempty(agentSet) 
                algorithm_indicator = 1;
                % PASS AGENT SET TO FORMATION CONTROLLER
                weights = [0.1;0.1;1;1];
                [v_boids] = obj.formationControl_boids(obj.targetWaypoint,agentSet,weights);
                % DEFINE VELOCITY REQUEST
                desiredVelocity = desiredSpeed*v_boids;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm      
            
            % GET THE EQUIVALENT HEADING ANGLE
            desiredHeadingVector = desiredVelocity/desiredSpeed;
            [dHeading,dPitch] = obj.GetVectorHeadingAngles([1;0;0],desiredHeadingVector); % Relative heading angles   
            omega_k_plus = [0;dPitch;-dHeading]/dt;         
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[desiredSpeed;0;0],omega_k_plus);
            obj.localState(1:6)  = obj.localState(1:6) + dt*dX;
            obj.localState(7:12) = dX;

            % ///////////// UPDATE OBJECT GLOBAL PROPERTIES ///////////////
            if obj.VIRTUAL.idleStatus
                obj.localState(7:12) = zeros(6,1);
            end
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,obj.localState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(TIME,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [obj.localState(7);obj.localState(4:6)];         % Record the control inputs
        end
    end
end