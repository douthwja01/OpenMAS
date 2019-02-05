%% POTENTIAL FIELD COLLISION AVOIDANCE AGENT (agent_potentialField2017.m) %%%%%%%%%%%%
% This agent uses a preliminary form of potential field theory to conduct
% its collision avoidance corrections.

% Author: James A. Douthwaite

classdef agent_potentialField2017 < agent
    
%% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties    
        % Performance Parameters
        maxAccelerations;
        minAccelerations;                   % Primative acceleration limits
        % Tuning parameters
        linearForceGradient;                % Linear avoidance force gradient (N/m)
        % gif toggle
        giffOn = 0;
    end
%%  CLASS METHODS
    methods 
        %% CONSTRUCTOR
        function obj = agent_potentialField2017(varargin)
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                           % Get super class 'agent'
            % AGENT SPECIFIC PARAMETERS
            obj.maxAccelerations = [20;20;20];
            obj.minAccelerations = [-20;-20;20];          % Primative acceleration limits
            obj.sampleFrequency = 19;                     % Virtual computation frequency
            obj.sensorRange = 50;                         % Virtual sensor range
            obj.linearForceGradient = 10;
            obj.mass = 1;                                 % Assign mass parameters (kg)
            % VIRTUAL DEFINITION
            obj.VIRTUAL.detectionRange = obj.sensorRange; % Assign the range attribute to the SIM VIRTUAL attribute
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
            
            % INPUT HANDLING
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            % CHECK FOR NEW INFORMATION UPDATE
            observationSet = varargin{1}; % The detected objects
            if isempty(observationSet)
                % NO INFORMATION AVAILABLE
                passiveStateUpdate = obj.stateDynamics(dt,[0;0;0],[0;0;0]); % Update the state vector
                obj = obj.updateGlobalProperties(dt,passiveStateUpdate);    % Update the objects global porperties
                return
            else
                % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
                [obj,obstacleSet,~] = obj.getAgentUpdate(observationSet);  
            end
            
            % PLOT THE LOCAL AVOIDANCE PROBLEM
            visualiseProblem = 0;
            if obj.objectID == 1 && visualiseProblem == 1
                overHandle = figure(100+obj.objectID);
                hold on; grid on;
                axis vis3d;
                view([-78 50]);
                xlabel('x_{m}');
                ylabel('y_{m}');
                zlabel('z_{m}');
                visualiseProblem = 1;
            end
            
            % 2. GET INSTANTANEOUS FORCES
            correctiveForce = zeros(3,1);
            for item = 1:length(obstacleSet)
                proximity = obstacleSet(item).state(1:3);
                [objectForce] = obj.linearForceModel(proximity,obj.linearForceGradient,visualiseProblem);
                % SUM THE FORCE CONTRIBUTIONS
                correctiveForce = correctiveForce + objectForce;
            end
            
            % 3. REVERT THE FORCE INTO A DESIRED ACCELERATION
            desiredAcceleration = correctiveForce/obj.mass;
            % CHECK AGAINST ACCELERATION CONSTRAINTS
            for i = 1:numel(desiredAcceleration)
                if desiredAcceleration(i) > obj.maxAccelerations(i)
                    desiredAcceleration(i) = obj.maxAccelerations(i);
                elseif desiredAcceleration(i) < obj.minAccelerations(i)
                    desiredAcceleration(i) = obj.minAccelerations(i);
                end
            end
            
            % 8. PLOT THE AVOIDANCE SITUATION
            if visualiseProblem && obj.objectID == 1
                % PLOT THE ESCAPE VELOCITY PROBLEM
                scatter3(desiredAcceleration(1)*dt,...
                         desiredAcceleration(2)*dt,...
                         desiredAcceleration(3)*dt,'r','filled');
                %  PLOT THE LOCAL AGENT VELOCITY
                q = quiver3(0,0,0,desiredAcceleration(1)*dt...
                                 ,desiredAcceleration(2)*dt...
                                 ,desiredAcceleration(3)*dt,'k');
                q.AutoScaleFactor = 1;
                view([-45 45 ]);
                drawnow;
                
                % ADD FIGURE TO GIFF
                fileName = sprintf('coneTransformation[%s]',obj.name);
                [obj] = obj.getAnimationFrame(overHandle,TIME,fileName);
                hold off;
                close;
            end
            
            % 9. UPDATE THE STATE VECTOR
            agentStateUpdate = obj.stateDynamics(dt,desiredAcceleration);
            % 10. UPDATE THE CLASS GLOBAL PROPERTIES
            obj = obj.updateGlobalProperties(dt,agentStateUpdate);
        end
    end
    
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        % GET LINEAR FORCE CONTRIBUTION FROM OBJECT 
        function [objectForce] = linearForceModel(obj,lambda,F_coeff,visualiseProblem)
            % This function derives the force contribution for a detected
            % obstacle using a linear force model.
            % INPUTS:
            % lambda           - The objects relative position [3x1] (m)
            % F_coeff          - The force coefficient for that object
            % OUTPUT:
            % objectForce      - The force contribution [3x1] (N)
            
            % INPUT HANDLING
             % IF NO Multiplier provided
            if ~exist('F_coeff','var')
                F_coeff = 1; % Assume a 1N-1m relationship
            end
            
%             display(lambda);
%             display(F_coeff);
            
            % Unitise the relative position
            mod_lambda = sqrt(sum(lambda.^2));
            unit_lambda = lambda/mod_lambda;
            
            % CALCULATE THE LINEAR FORCE CONTRIBUTION 
            objectForce = (F_coeff/mod_lambda).*-unit_lambda;              % Force acts towards the agent
                       
            display(objectForce);
            % PLOT THE PROBLEM IF REQUIRED
            if exist('visualiseProblem','var') && visualiseProblem == 1
                proxSeries = 0:lambda/10:lambda(1);
                forceSeries = 0:objectForce/10:objectForce;
                plot(forceSeries,proxSeries);
%                 %  PLOT THE LOCAL AGENT VELOCITY
%                 q = quiver3(0,0,0,desiredAcceleration(1)*dt...
%                                  ,desiredAcceleration(2)*dt...
%                                  ,desiredAcceleration(3)*dt,'k');
%                 q.AutoScaleFactor = 1;
                
            end
            

        end
    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Static)        

    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]