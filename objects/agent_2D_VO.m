%% 2D GEOMETRIC COLLISION AVOIDANCE AGENT (agent_VO_2D %%%%%%%%%%%%%%%%%%%%


% Author: James A. Douthwaite

classdef agent_2D_VO < agent_2D & agent_VO
    properties
    end
    
    methods
        % CONSTRUCTION METHOD
        function obj = agent_VO_2D(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        % AGENT MAIN CYCLE
        function obj = processTimeCycle(obj,TIME,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            observationSet = varargin{1};                                  % The detected objects
            % The agent_2D provides a 2D interface for the 3D object world,
            % obstacle states are instead of the form [x y u v]. 
            [obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(observationSet); 

            % GET THE TARGET WAYPOINT (DESIRED) HEADING
            if ~isempty(obj.targetWaypoint)  
               targetHeading = obj.targetWaypoint.state(1:2)/norm(obj.targetWaypoint.state(1:2));
            else
               targetHeading = [1;0];
            end
            targetSpeed = 2;

            % ////////////////// OBSTACLE AVOIDANCE ///////////////////////
            avoidanceSet = [obstacleSet,agentSet];
            algorithm_start = tic; algorithm_indicator = 0;  avoidanceEnabled = 1;  
            if ~isempty(avoidanceSet) && avoidanceEnabled
                algorithm_indicator = 1;
                % GET THE UPDATED DESIRED VELOCITY
                [desiredVelocity] = obj.getAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);
            newSpeed = (norm(obj.localState(7:8)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = d_heading;             
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            newState = obj.stateDynamics_velocities(dt,newSpeed,newHeading);
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:2,TIME.currentStep) = [newSpeed;newState(12)];      % Record the control inputs 
            obj.DATA.inputNames = {'Speed (m/s)','Yaw (rad)'};
        end
        
        % CALCULATE THE NECESSARY 2D AVOIDANCE VELOCITY
        function [desiredVelocity] = obj.getAvoidanceCorrection(dt,desiredVelocity,avoidanceSet,visualiseProblem);
            
        end
        
    end
end
% AGENT STATE VECTOR [x;y;v;u;psi;the;phi;p;q;r]