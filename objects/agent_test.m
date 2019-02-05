%% TEST AGENT 
% This agent is used to examine the simulation feedback.

% Author: James A. Douthwaite

classdef agent_test < agent 
    
    properties
    end
    %  CLASS METHODS
    methods 
        % CONSTRUCTOR
        function obj = agent_test(varargin)
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Get super class 'agent'

%             obj.VIRTUAL.detectionRadius = 10;
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % AGENT MAIN CYCLE 
        function [obj] = main(obj,ENV,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('Object TIME packet is invalid.');
            end           

            % DEFAULT BEHAVIOUR
            desiredVelocity = [1;0;0]*obj.nominalSpeed;
            
            % UPDATE THE AGENT WITH THE NEW INFORMATION
            [obj,obstacleSet,agentSet,waypointSet] = obj.getAgentUpdate(varargin{1});  
                            
            % PLOTTING THE COMPLETE OBJECT SURROUNDINGS
            if obj.objectID == 1
                figureHandle = figure(1);
                [figureHandle] = obj.getObjectScene([obstacleSet;agentSet;waypointSet],figureHandle);   % Plot the objects
                [figureHandle] = obj.getAnimationFrame(ENV,figureHandle,'objectScene.gif');             % Store the annimations
            end
            
%             % GET THE COMPLEX OBSTACLE REPRESENTATION
%             for i = 1:numel(obstacleSet)
%                 [VO_complex] = obj.getComplexObstacles(zeros(3,1),...
%                     obj.localState(7:9),...
%                     obj.VIRTUAL.radius,...
%                     obstacleSet.position,...
%                     obstacleSet.velocity,...
%                     obstacleSet.geometry);
%             end
            
            
            
            
            
            
            
            % REQUEST A SEQUENCE OF ROTATIONAL RATES
            [omega] = obj.proceduralHeadingRate(ENV,3);
%             omega = zeros(3,1);
            % REQUEST A CONSTANT SPEED
            velocity = [norm(desiredVelocity);0;0];
            
            % ///////////////////// STATE DYNAMICS ////////////////////////
            [newState] = obj.updateLocalState(ENV,obj.localState,velocity,omega);
            
            % //////////////// UPDATE GLOBAL PROPERTIES ///////////////////
            [obj] = obj.updateGlobalProperties_TEST(dt,newState);
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            obj.DATA.inputNames = {'Phi (rad)','Pitch (rad)','Psi (rad)'};
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
        % INITIALISE THE STATE VECTOR AS [x y z phi theta psi]
        function [obj] = initialise_localState(obj,localXYZvelocity,localXYZrotations)
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
        % GLOBAL UPDATE - EULER 6DOF(3DOF) [x y z phi theta psi],[x y psi]
        function [obj] = updateGlobalProperties_TEST(obj,dt,eulerState)
            % This function computes the position, velocity and 
            % quaternion in the global ENU frame from a 6DOF state vector 
            % with euler rotations.
            
            % This function is intended for state vectors of the form:
            % - [x y z phi theta psi]
            % - [x y phi]
            
            % DEFINE ROTATION INDICES FOR 2D AND 3D ROTATING SYSTEMS
            if numel(eulerState) == 6
                positionIndices = 1:3;
                eulerIndices = 4:6;                     % The rotations about X,Y,Z
            elseif numel(eulerState) == 3
                positionIndices = 1:2;
                eulerIndices = 3;                       % The rotations about Z
            else
                error('State notation not recognised');
            end
            
            % EQUIVALENT RATES
            velocity_k_plus   = (eulerState(positionIndices) - obj.VIRTUAL.priorState(positionIndices))/dt;
            eulerRates_k_plus = (eulerState(eulerIndices) - obj.VIRTUAL.priorState(eulerIndices))/dt;  
            
            % ///// IF ALL WAYPOINTS ARE ACHEIVED; FREEZE THE AGENT ///////
            if isprop(obj,'targetWaypoint')
                if isempty(obj.targetWaypoint) && ~isempty(obj.achievedWaypoints)
                    obj.VIRTUAL.idleStatus = logical(true);
                    velocity_k_plus   = zeros(numel(positionIndices),1);   % Freeze the agent
                    eulerRates_k_plus = zeros(numel(eulerIndices),1);      % Freeze the agent
                end   
            end

            % EULER ROTATIONS -> GLOBAL ROTATIONS
            if numel(eulerRates_k_plus) ~= 3
                localAxisRates = [0;0;1]*eulerRates_k_plus; 
                velocity_k_plus = [velocity_k_plus;0];
            else
                % CALCULATE THE GLOBAL AXIS RATES (FROM 3DOF)
                localAxisRates = eulerRates_k_plus;                        % Rotations rates are negative
            end
            
            % /////////// DISPLAY OBJECT ID FOR UPDATE/CLEARITY ///////////
%             fprintf('Updating %s [ID:%s]\n',obj.name,num2str(obj.objectID));
                        
            % MAP LOCAL RATES TO THE NEW QUATERNION POSE (maps a local vector to global)
            % quaternion_k_plus & obj.VIRTUAL.quaternion must define the 
            % rotation from the global axes to the new body pose. 

%             % PREVIOUS QUATERNION              % [ IF THE QUATERNION DEFINES B-G ]   
%             quaternion_k = obj.VIRTUAL.quaternion;            
%             % THE LOCAL AXIS RATES              
%             omega = eulerRates_k_plus;
%             % UPDATE THE QUATERNION POSE
%             quaternion_k_plus = OMAS_geometry.integrateQuaternion(quaternion_k,omega,dt)
%             % REDEFINE THE ROTATION MATRIX
%             R_k_plus = quat2rotm(quaternion_k_plus')
            
%             display('test values');

            % MAP THE LOCAL RATES TO GLOBAL RATES AND INTEGRATE QUATERNION
            % PREVIOUS QUATERNION                % [ IF THE QUATERNION DEFINES GB ]
            quaternion_k = obj.VIRTUAL.quaternion;   
            % PREVIOUS ROTATION-MATRIX
            R_k_plus = quat2rotm(quaternion_k');
            % THE GLOBAL AXIS RATES       
            omega = R_k_plus'*localAxisRates;
            % UPDATE THE QUATERNION POSE
            quaternion_k_plus = OMAS_geometry.integrateQuaternion(quaternion_k,omega,dt);
            % REDEFINE THE ROTATION MATRIX
            R_k_plus = quat2rotm(quaternion_k_plus');
            
%             % MAP EULER ANGLE STATE TO NEW ROTATION MATRIX
%             angles = eulerState(4:6);
%             R_k_plus = eul2rotm(angles','XYZ');
%             quaternion_k_plus = rotm2quat(R_k_plus)'
%             quaternion_k_plus = eul2quat(angles','XYZ')';
%             R_k_plus = quat2rotm(quaternion_k_plus');

            % MAP THE LOCAL VELOCITY TO THE GLOBAL AXES
            globalVelocity_k_plus = R_k_plus'*velocity_k_plus;
            globalPosition_k_plus = obj.VIRTUAL.globalPosition + dt*obj.VIRTUAL.globalVelocity;
            % ///////////////// REASSIGN K+1 PARAMETERS ///////////////////
            [obj] = obj.updateGlobalProperties_direct(globalPosition_k_plus,... % Global position at k plius
                                                      globalVelocity_k_plus,... % Global velocity at k plus
                                                      quaternion_k_plus,...     % Quaternion at k plus
                                                      eulerState);              % The new state for reference
        end
    end
end