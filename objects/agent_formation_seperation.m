%% FORMATION CONTROL AGENT (SEPERATION) (agent_formation_seperation.m) %%
% This agent is designed to implement the seperation-based control law
% used in our current coordination studies. The agent has the ability to
% design and inact a control vector that drives the agents towards a
% formation defined by an adjacency matrix.

% Author: James A. Douthwaite & Dr Shiyu Zhao

classdef agent_formation_seperation < agent
    % INITIALISE THE PROPERTIES OF THE FORMATION CONTROLLER
    properties
        % FORMATION PARAMETERS
        adjacencyMatrix;
    end
    %  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_formation_seperation(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % GET THE VELOCITCY OBSTACLE (VO) AGENT TOOLS
            obj@agent(varargin);                                       % Get the supercalss
            
            % KINEMATIC LIMITS
            obj.linearVelocityLimits  = 2*ones(3,1);
%             obj.angularVelocityLimits = 10*ones(3,1);
            obj.priorError = zeros(4,1);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);
        end
        
        %% AGENT MAIN CYCLE
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % PLOT AGENT FIGURE
            visualiseProblem = 0;
            visualiseAgent = 1;
            if obj.objectID == visualiseAgent && visualiseProblem == 1
                overHandle = figure(100 + obj.objectID);
                hold on; grid on;
                xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
            end 
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            observationSet = varargin{1};                                  % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,agentSet,~] = obj.getAgentUpdate(observationSet);
            desiredSpeed = 0;
            desiredVelocity = [1;0;0]*desiredSpeed;                        
            
            % ////////////////// FORMATION CONTROLLER /////////////////////
            % We wish to conduct formation control with the other agents
            % in the simulation only.
            algorithm_start = tic; algorithm_indicator = 0;  formationEnabled = 1;  
            if ~isempty(agentSet) && formationEnabled
                algorithm_indicator = 1;
                % PASS AGENT SET TO FORMATION CONTROLLER
                [fi,V] = obj.formationControl_seperation(agentSet);        % Get the force vector
                % NORMALISE THE FORCE VECTOR
                norm_fi = norm(fi);
                unit_fi = fi/norm_fi; % < -------------------- TO BE CHANGED
                constrainedNorm = obj.boundValue(norm_fi,-obj.linearVelocityLimits(1),obj.linearVelocityLimits(1));
                % DEFINE VELOCITY REQUEST
                desiredVelocity = constrainedNorm*unit_fi;
            end
            algorithm_dt = toc(algorithm_start);                           % Stop timing the algorithm                      
            
            % /////////////// PARSE CONTROLLER INPUTS /////////////////////
            targetSpeed = norm(desiredVelocity);
            targetHeading = desiredVelocity/targetSpeed;                   % Get the control inputs
            
            % ////////////////// AGENT CONTROLLER /////////////////////////
            [d_heading,d_speed,obj] = obj.simpleController(targetHeading,targetSpeed);              
            newSpeed = (norm(obj.localState(7:9)) + d_speed);              % Define the velocity vector from the new speed
            newHeading = zeros(3,1) + d_heading;                                        % Heading change relative to straight ahead
            
            % //////////////// APPLY KINEMATIC LIMITS /////////////////////
            [newHeading,newSpeed] = obj.kinematicContraints(dt,newHeading,newSpeed);
            newVelocity = [1;0;0]*newSpeed;                                % Convert to local velocity vector
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            [newState] = obj.stateDynamics_velocities(dt,newVelocity,newHeading);
            % IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT
            if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                newState = obj.freezeAgent();
            end
            
            % /////////// UPDATE THE CLASS GLOBAL PROPERTIES //////////////
            obj = obj.updateGlobalProperties(dt,newState);
            
            % ////////////////// UPDATE AGENT STATE ///////////////////////
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:length(obj.priorError),TIME.currentStep) = [newSpeed;newState(4:6)];         % Record the control inputs 
            obj.DATA.inputNames = {'Speed (m/s)','Roll (rad)','Pitch (rad)','Yaw (rad)'};
            obj.DATA.lypanov(TIME.currentStep) = V;                               % Record the lypanov value
        end

    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Access = public)
        % FORMATION CONTROLLER
        function [fi,V] = formationControl_seperation(obj,observedAgents)
            % This function calculates the formation control vector to
            % to bring about the desired object seperation, given the 
            % agents current knowledge of the surrounding agents.
            
            % Author: Dr Shiyu Zhao
            
            % Under this implementation, the adjacency matrix is indexed by
            % agent objectID's. This means agent 2 will only calculate the
            % contributions from the observed agents 1 and 2.
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedAgents);
            pi = obj.localState(1:3);                                    % Other objects are observed relatively
            
            % CONFIRM ADJACENCY MATRIX
            if ~isprop(obj,'adjacencyMatrix')
                error('Agent is missing required adjacency matrix');
            end
            
            fi = zeros(3,1);
            V = 0;
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            for j = 1:agentNumber
                % SECOND AGENT INDEX
                agentID = observedAgents(j).objectID;
                pj = pi + observedAgents(j).state(1:3);  
                % SEPERATION BASED CONFIRMATION
                ell_ij = obj.adjacencyMatrix(obj.objectID,agentID);
                % SUM THE CONTRIBUTION
                fi = fi + (norm(pi - pj)^2 - ell_ij^2)*(pj - pi);  
                % GET THE LYAPANOV VALUE
                V = V + (norm(pi-pj)^2-ell_ij^2)^2;
            end
        end
        % POSITION BASED FORMATION CONTROLLER
        function [fi,V] = formationControl_position(obj,observedAgents)
            % This function calculates the formation control vector to
            % to bring about the desired object position, given the 
            % agents current knowledge of the surrounding agents.
            
            % Author: Dr Shiyu Zhao
            
            % Under this implementation, the adjacency matrix is indexed by
            % agent objectID's. This means agent 2 will only calculate the
            % contributions from the observed agents 1 and 2.
            
            % CONSTANT PARAMETERS
            agentNumber = numel(observedAgents);
            pi = obj.localState(1:3);                                    % Other objects are observed relatively
                     
            fi = zeros(3,1);
            V = 0;            
            % MOVE THROUGH ALL AGENTS AN CALCULATE THEIR FORCE CONTRIBUTION
            % We are i
            
            for j = 1:agentNumber
                % SECOND AGENT INDEX
                agentID = observedAgents(j).objectID;
                pj = pi + observedAgents(j).state(1:3);  
                pj = p_all(:,j);                                                   % Position of the second agents from the global set
                
                pi_star = p_star_all(:,i);
                pj_star = p_star_all(:,j);
                
                % SEPERATION BASED CONFIRMATION
                ell_ij = obj.adjacencyMatrix(obj.objectID,agentID);
                % SUM THE CONTRIBUTION
                fi = fi + (pj - pi - (pj_star - pi_star));
                % GET THE LYAPANOV VALUE
                V = V + (norm(pi-pj)^2-ell_ij^2)^2;
            end 
            

        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]