%% THE 2D AGENT BASE CLASS (agent.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to the 'agent' class, this object class is intended provide
% support for 2D simulation. The functions provided here 

% Author: James A. Douthwaite

classdef agent_2D_test < agent_2D
%%% AGENT(2D) BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties

    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = agent_2D_test(varargin)
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
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent_2D(varargin); 
            
            obj.localState = zeros(4,1);
            
            % AGENT DEFAULT KINEMATIC LIMITATIONS
%             obj.linearVelocityLimits = [inf;inf];
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % ///////////////// AGENT MAIN CYCLE //////////////////////////////
        function obj = main(obj,TIME,varargin)
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
            observationSet = varargin{1}; % The detected objects
            [obj,obstacleSet,agentSet,waypointSet] = obj.GetAgentUpdate(dt,observationSet);

            % GET WAYPOINT (TARGET) HEADING VECTOR \\\\\\\\\\\\\\\\\\\\\\\\
            targetHeading = obj.GetTargetHeading();
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            targetSpeed = 2;
            desiredVelocity = targetHeading*targetSpeed;
            
            % PASS THE REQUEST TO THE CONTROLLER
            [obj] = obj.controller(dt,desiredVelocity);
                   
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            obj.DATA.inputNames = {'dx (m/s)','dy (m/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),TIME.currentStep) = [obj.localState(3:4,1)];         % Record the control inputs 
        end
        % /////////////////////////////////////////////////////////////////

        % UPDATE GLOBAL PROPERTIES
        function [obj] = updateGlobalProperties(obj,dt,localState_update)
           % This function computes the new global parameters for the given
           % agent based on its new state.
           
            % CONSTANT PARAMETERS
            velocityIndices = 3:4;
            
            % DEFINE UPDATE PARAMETERS
            globalPosition_k = obj.VIRTUAL.globalPosition;                 % 3D although the 2D
            globalVelocity_k = obj.VIRTUAL.globalVelocity;
            quaternion_k     = obj.VIRTUAL.quaternion;
            
            localState_ENU_k = obj.localState; %<---- must be redefined
            localState_ENU_k_plus = localState_update; 
            
%             localState_FLU_k = zeros(4,1);
%             localState_FLU_k(1) =  localState_FRD_k(1);
%             localState_FLU_k(2) = -localState_FRD_k(2);
%             localState_FLU_k(3) =  localState_FRD_k(3);
%             localState_FLU_k(4) = -localState_FRD_k(4)
            
%             localState_FLU_k_plus = zeros(4,1);
%             localState_FLU_k_plus(1) =  localState_FRD_k_plus(1);
%             localState_FLU_k_plus(2) = -localState_FRD_k_plus(2);
%             localState_FLU_k_plus(3) =  localState_FRD_k_plus(3);
%             localState_FLU_k_plus(4) = -localState_FRD_k_plus(4)
            
            localState_XYZ_k = localState_ENU_k;
            localState_XYZ_k(2) = localState_XYZ_k(2);
            localState_XYZ_k(4) = localState_XYZ_k(4);
            
            localState_XYZ_k_plus = localState_ENU_k_plus;
            localState_XYZ_k_plus(2) = localState_XYZ_k_plus(2);
            localState_XYZ_k_plus(4) = localState_XYZ_k_plus(4);
                       
            temp_k = [localState_XYZ_k(3:4,1);0];
            temp_k_plus = [localState_XYZ_k_plus(3:4,1);0];
            
            % GET THE NEW ROTATION QUATERNION 
            [dq] = OMAS_geometry.vectorsToQuaternion(temp_k_plus,temp_k);
            quaternion_k_plus = OMAS_geometry.qMultiply(dq,quaternion_k);
            
            
            % NEW ROTATION MATRIX
            [R_k_plus] = OMAS_geometry.quaternionToRotationMatrix(quaternion_k_plus);
            % ////////// UPDATE GLOBAL POSE ///////////////////////////////
            globalVelocity_k_plus = [localState_XYZ_k_plus(velocityIndices);0];
            globalPosition_k_plus = globalPosition_k + dt*globalVelocity_k;
            
            % ////////// REASSIGN K+1 PARAMETERS //////////////////////////
            obj.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            obj.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            obj.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            obj.VIRTUAL.rotationMatrix = R_k_plus;
            obj.localState = localState_update;  
        end
    end
    %% STATIC AGENT TOOLS
    methods (Static)
       
    end  
end

