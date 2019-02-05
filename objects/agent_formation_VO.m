%% VELOCITY OBSTACLE/FORMATION CONTROL AGENT (agent_formation_VO.m) %%%%%%%
% This agent is a hybrid of the formation control element agent class and
% the Velocity Obstacle (VO) class. 

% Author: James A. Douthwaite

classdef agent_formation_VO < agent_VO & agent_formation
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_formation_VO(varargin)

            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_VO(varargin);
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);   
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
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
                [vi,L] = obj.formationControl_distance(agentSet);        % Get the force vector
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
                [desiredHeadingVector,desiredSpeed] = obj.getAvoidanceCorrection(desiredVelocity,avoidanceSet,visualiseProblem);
                desiredVelocity = desiredHeadingVector*desiredSpeed;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm 
            % APPLY SPEED CONSTRAINT
            if norm(desiredVelocity) > obj.maxSpeed
                desiredSpeed = obj.maxSpeed;
                desiredVelocity = (desiredVelocity/norm(desiredVelocity))*desiredSpeed;
            end   
            
            % //////////////////// AGENT CONTROL //////////////////////////
            % APPLY SPEED CONSTRAINT
            desiredSpeed = norm(desiredVelocity);
            if desiredSpeed > obj.maxSpeed
                desiredHeadingVector = desiredVelocity/norm(desiredVelocity);
                desiredVelocity = desiredHeadingVector*obj.maxSpeed;
            elseif desiredSpeed == 0
                desiredHeadingVector = [1;0;0];
            end  
            
            % GET THE EQUIVALENT HEADING ANGLE
            [dHeading,dPitch] = obj.getVectorHeadingAngles([1;0;0],desiredHeadingVector); % Relative heading angles   
            omega_k_plus = [0;dPitch;-dHeading]/dt;           
            
            % /////////// SIMPLE DYNAMICS + PURE TRANSLATION //////////////
            [dX] = obj.dynamics_simple(obj.localState(1:6),[desiredSpeed;0;0],omega_k_plus);
            obj.localState(1:6)  = obj.localState(1:6) + dt*dX;
            obj.localState(7:12) = dX;
              
            % ////// GLOBAL UPDATE FOR STATE WITH RETAINED VELOCITES //////
            [obj] = obj.updateGlobalProperties_3DVelocities(dt,obj.localState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj = obj.writeAgentData(TIME,algorithm_indicator,algorithm_dt);
            obj.DATA.inputNames = {'Vx (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [obj.localState(7);obj.localState(4:6)];         % Record the control inputs
            obj.DATA.lypanov(TIME.currentStep) = L;                               % Record the lypanov value
        end
    end
    methods
        % INITIALISE THE [x;x_dot]' 3D STATE VECTOR (FOR ALL 'VO' SUBCLASSES
        function [obj] = initialise_localState(obj,localXYZVelocity,localXYZrotations)
            % The state initialiser must be called 'initialise_localState'
            % and instead calls the 'initialise_3DVelocities' function in
            % this case. 
            
            % BUILD THE STATE VECTOR FOR A 3D SYSTEM WITH CONCATINATED VELOCITIES
            [obj] = obj.initialise_3DVelocities(localXYZVelocity,localXYZrotations);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]