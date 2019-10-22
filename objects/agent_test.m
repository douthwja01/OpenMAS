%% TEST AGENT 
% This agent is used to examine the simulation feedback.

% Author: James A. Douthwaite

classdef agent_test < agent 
    
    properties
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function obj = agent_test(varargin)

            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'
            
            % Assign defaults
            obj.localState = zeros(6,1);
            [obj] = obj.SetBufferSize(5);
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end  
        % Setup
        function [obj] = setup(obj,localXYZVelocity,localXYZrotations)
            % This function is called in order to build the initial state
            % vector for the generic agent class 'objectDefinition'.
            % INPUTS:
            
            % ASSUMPTIONS:
            % - The designed state is described in the ENU axes.
            % - The state is of the form:
            % - [x y z phi theta psi]
            
            % BUILD STATE VECTOR
            obj.localState = zeros(6,1);
            obj.localState(4:6) = localXYZrotations;  % The euler rotations are rotations about the local X,Y,Z respectively
            % RETAIN THE PRIOR STATE FOR REFERENCE
            obj.VIRTUAL.priorState = obj.localState;
        end
        % Main
        function [obj] = main(obj,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
                 

%             % PLOT AGENT FIGURE
%             visualiseProblem = 1;
%             visualiseAgent = 1;
%             if obj.objectID == visualiseAgent && visualiseProblem == 1
%                 overHandle = figure('name','testFigure');
%                 ax = axes(overHandle);
%                 hold on; grid on;
%                 axis equal;
%                 xlabel('x_{m}'); ylabel('y_{m}'); zlabel('z_{m}');
%             end 
            
            
            % DEFAULT BEHAVIOUR
            desiredVelocity = [1;0;0]*obj.nominalSpeed;
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(ENV,varargin{1});
            % /////////////////////////////////////////////////////////////
            
            % REQUEST A SEQUENCE OF ROTATIONAL RATES
%             [omega] = obj.proceduralHeadingRate(ENV,3);
            omega = zeros(3,1);
            % REQUEST A CONSTANT SPEED
            velocity = [norm(desiredVelocity);0;0];
            
            % ///////////////////// STATE DYNAMICS ////////////////////////
            [newState] = obj.updateLocalState(ENV,obj.localState,velocity,omega);
            
            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            [obj] = obj.updateGlobalProperties(ENV.dt,newState);
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            obj.DATA.inputNames = {'$\phi (rad)$','$\theta (rad)$','$\psi (rad)$'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = newState(4:6);         % Record the control inputs
        end       
    end

    methods
        % DECOMPOSE THE OBSTACLE GEOMETRIES INTO A VO SET
        function [VO] = getComplexObstacles(obj,p_a,v_a,r_a,p_b,v_b,geometry,tau)
            % This function computes the VO set from a list of complex
            % obstacles defined by their .GEOMETRY parameter.
            
            % NOTES:
            % - The .geometry parameter is the geometry of the object pre-
            %   positioned relative to the agent, rotated and scaled.
            % - VO structure:
            %   VO = struct('apex',v_b,...
            %     'axisUnit',axisUnit,...
            %     'axisLength',mod_lambda_ab,...
            %     'openAngle',2*halfAlpha,...
            %     'leadingEdgeUnit',unit_leadingTangent,...
            %     'trailingEdgeUnit',unit_trailingTangent,...
            %     'isVaLeading',isVaLeading,...
            %     'isVaInsideCone',0,...
            %     'truncationTau',tau,...
            %     'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
            %     'truncationCircleRadius',(r_a + r_b)/tau);
            
            norm_lateralPositionProjection = norm([p_b(1:2);0]);
            unit_lateralPositionProjection = [p_b(1:2);0]/norm_lateralPositionProjection;   % XY PLANAR PROJECTION OF RELATIVE POSITION
            
            a_min = 0; a_max = 0;
            % Get the points with the largest perpendicular projection
            for v = 1:size(geometry.vertices,1)
                % GET THE XY VERTEX PROJECTION
                lateralVertexProjection = [geometry.vertices(v,1:2)';0];          % Projection of the vertex in the XY plane
                unit_lateralVertexProjection = lateralVertexProjection/norm(lateralVertexProjection);
                % GET THE ANGULAR ARGUEMENT WITH THE POSITION VECTOR
                crossProduct = cross(unit_lateralVertexProjection,unit_lateralPositionProjection); % Planar normal
                dotProduct   = dot(unit_lateralVertexProjection,unit_lateralPositionProjection);
                % THE SIGNED ANGLE OF THE VERTEX RELATIVE TO THE
                % CENTROID VECTOR
                signedAngularProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
                % RETAIN THE LARGEST HEADING CHANGE AS THE VO BOUNDARY
                if signedAngularProjection > a_max
                    a_max = signedAngularProjection;
                    unit_trailingTangent = unit_lateralVertexProjection; % Assuming that clockwise +ve
%                        maximalVertexProjection = [obstacleGeometry.vertices(v,1:2)';0];
                elseif signedAngularProjection < a_min
                    a_min = signedAngularProjection;
                    unit_leadingTangent = unit_lateralVertexProjection;  % Assuming that anti-clockwise -ve
%                        minimalVertexProjection = [obstacleGeometry.vertices(v,1:2)';0];
                end
            end
            
            % CHECK IF Va belongs to the projection cone
            crossProduct = cross(v_a,unit_lateralPositionProjection);      % Planar normal
            dotProduct = dot(v_a/norm(v_a),unit_lateralPositionProjection);% Angle between two
            signedVelocityProjection = atan2(dot(crossProduct,[0;0;1]),dotProduct);
            % DETERMINE IF THE CURRENT VELOCITY BELONGS TO THE CURRENT VO
            isVaInsideCone = 0;
            if signedVelocityProjection < a_max && signedVelocityProjection > a_min
                isVaInsideCone = 1;
            end
            
            % THE MAXIMAL VERTEX PROJECTIONS
            equivalentOpenAngle = abs(a_min) + abs(a_max);
            isVaLeading = 1;                                               % Va always leading (obstacle is static)

            effectiveRadii = norm_lateralPositionProjection*sin(equivalentOpenAngle);       % <<< TO BE DEFINED         %(r_a + r_b)
            % DEFINE THE VO STRUCTURE
            VO = struct('apex',v_b,...
                'axisUnit',unit_lateralPositionProjection,...
                'axisLength',norm_lateralPositionProjection,...
                'openAngle',equivalentOpenAngle,...
                'leadingEdgeUnit',unit_leadingTangent,...
                'trailingEdgeUnit',unit_trailingTangent,...
                'isVaLeading',isVaLeading,...
                'isVaInsideCone',isVaInsideCone,...
                'truncationTau',tau,...
                'truncationCircleCenter',(p_b - p_a)/tau + v_b,...
                'truncationCircleRadius',(effectiveRadii)/tau);            % Double the radii of A clear
        end
    end
    
    methods (Static)
        % PROCEDURAL MOVEMENT FUNCTION
        function [omega_k] = proceduralHeadingRate(TIME,cycleTime)
            % This function requests a sequence of rotational rates to test
            % the mapping of states-euler angles to the global state
            % rotations.
            
            omega_k = [0;0;0]; % -ve feedback
            if TIME.currentTime >= 0 && TIME.currentTime < cycleTime
                omega_k(1) = 0.5;
                fprintf('rolling...\n');
            end
            if TIME.currentTime >= cycleTime && TIME.currentTime < cycleTime*2
                omega_k(2) = 0.5;
                fprintf('pitching...\n');
            end
            if TIME.currentTime >= cycleTime*2 && TIME.currentTime < cycleTime*3
                omega_k(3) = 0.5;
                fprintf('yawing...\n');
            end
        end
    end
    %% ///////////////////////// TEST FUNCTIONS ///////////////////////////
    methods
        % GET THE STATE UPDATE (USING ODE45)
        function [X] = updateLocalState(obj,TIME,X,velocity,omega)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                return
            else
                tspan = [TIME.currentTime TIME.timeVector(TIME.currentStep + 1)];
                opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
                [~,Xset] = ode45(@(t,X) obj.dynamics_simple(X,velocity,omega), tspan, X, opts);
                X = Xset(end,:)';
            end
        end
    end
end