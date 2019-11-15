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
    
    %% ////////////////////////// MAIN METHODS ////////////////////////////
    methods 
        % Constructor
        function [this] = agent_2D_test(varargin)
            % This function is to construct the agent object using the
            % object defintions held in the 'objectDefinition' base class.
            % INPUTS:
            % namestr - Assigned name string
            % OUTPUTS:
            % obj     - The constructed object
            
            % Call the super class
            this@agent_2D(varargin); 
            
            % Assign defaults
            this.localState = zeros(4,1);
            
            % //////////////// Check for user overrides ///////////////////            
            % - It is assumed that overrides to the properties are provided
            %   via the varargin structure.
            [this] = this.ApplyUserOverrides(varargin); 
            % /////////////////////////////////////////////////////////////
        end
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % TIME     - The TIME simulations structure
            % varargin - Cell array of inputs
            % OUTPUTS:
            % obj      - The updated project
                        
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            observationSet = varargin{1}; % The detected objects
            [this,obstacleSet,agentSet,waypointSet] = this.GetAgentUpdate(ENV,observationSet);

            % GET WAYPOINT (TARGET) HEADING VECTOR \\\\\\\\\\\\\\\\\\\\\\\\
            targetHeading = this.GetTargetHeading();
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                             %
            %                         DO NOTHING                          %
            %                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ASSUME THERE ARE NO INERTIAL ACCELERATIONS
            targetSpeed = 2;
            desiredVelocity = targetHeading*targetSpeed;
            
            % PASS THE REQUEST TO THE CONTROLLER
            [this] = this.Controller(ENV.dt,desiredVelocity);
                   
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this.DATA.inputNames = {'dx (m/s)','dy (m/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = [this.localState(3:4,1)];         % Record the control inputs 
        end
    end
    methods
        % UPDATE GLOBAL PROPERTIES
        function [this] = UpdateGlobalProperties_TEST(this,dt,localState_update)
           % This function computes the new global parameters for the given
           % agent based on its new state.
           
            % CONSTANT PARAMETERS
            velocityIndices = 3:4;
            
            % DEFINE UPDATE PARAMETERS
            globalPosition_k = this.VIRTUAL.globalPosition;                 % 3D although the 2D
            globalVelocity_k = this.VIRTUAL.globalVelocity;
            quaternion_k     = this.VIRTUAL.quaternion;
            
            localState_ENU_k = this.localState; %<---- must be redefined
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
            this.VIRTUAL.globalVelocity = globalVelocity_k_plus;            % Reassign the global velocity
            this.VIRTUAL.globalPosition = globalPosition_k_plus;            % Reassign the global position
            this.VIRTUAL.quaternion = quaternion_k_plus;                    % Reassign the quaternion
            this.VIRTUAL.rotationMatrix = R_k_plus;
            this.localState = localState_update;  
        end
    end
end

