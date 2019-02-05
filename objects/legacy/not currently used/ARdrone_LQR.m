%% ARDRONE TEST

% Author: James A. Douthwaite

classdef ARdrone_LQR < ARdrone
    % LQR SPECIFIC PARAMETERS
    properties
        K_lqr;
    end

    methods 
        %% CONSTRUCTOR METHOD
        function obj = ARdrone_LQR(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 

            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@ARdrone(varargin);                                         % Create the super class 'agent'                  
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(varargin);            
            
            % DEFINE THE LQR PROPERTIES
            obj.K_lqr = obj.getLQRFeedback(obj.SYS.A,obj.SYS.B);           % Get the descrete LQR feedback   
        end
        %% AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % CHECK FOR NEW INFORMATION UPDATE \\\\\\\\\\\\\\\\\\\\\\\\\\\\
            observationSet = varargin{1}; % The detected objects
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,~] = obj.getAgentUpdate(observationSet);
            
            % INPUT STEP
            stateReference = zeros(12,1);
            if TIME.currentTime > 2
                stateReference = [2;0;-1;0;0;0;0;0;0;0;0;0];
            end
            %% ////// CONTROLLER //////////////////////////////////////////
            algorithm_start = tic; algorithm_indicator = 1;
            
            % CALCULATE THE STATE ERROR
            stateError = stateReference - obj.localState;
            % GET THE PERTURBATION INPUTS
            inputs_dU = obj.K_lqr*stateError; 
            
            % GET THE INPUTS REQUIRED TO HOVER
            [inputs_ss] = obj.inverseKinematics([-9.81;0;0;0]);    % NOMINAL INPUTS
            % SUM FOR THE ABSOLUTE INPUTS
            inputs_abs = inputs_ss + inputs_dU;            
            % END TIMER
            algorithm_dt = toc(algorithm_start);
            
            % STATE UPDATE FUNCTION INHERITED FROM 'ARdrone' //////////////
            [newState] = obj.stateDynamics(obj.localState,dt,inputs_abs);           
            % GENERIC STATE UPDATE FUNCTION
            obj = obj.updateGlobalProperties(dt,newState);
            % \\\\\\\\\\\\\\\ RECORD THE AGENT-SIDE DATA \\\\\\\\\\\\\\\\\\
            obj.DATA.algorithm_indicator(TIME.currentStep) = algorithm_indicator; % Record when the algorithm is ran
            obj.DATA.algorithm_dt(TIME.currentStep) = algorithm_dt;               % Record the computation time
            obj.DATA.inputs(1:4,TIME.currentStep) = inputs_abs;                   % Record the control inputs 
        end
        
    end
    
    methods (Static)
        % GET LINEAR-QUADRATIC-REGULATOR(LQR) CONTROLLER
        function [K_lqr] = getLQRFeedback(A,B)
            % This function designs the LQR controller used to provide
            % error feedback.
            % INPUTS:
            % A - The descrete plant matrix (dt_sim)
            % B - The descrete input matrix (dt_sim)
            
            % SIMULINK CONTROLLER AND PARAMETERS
%             Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1E1;
%             R = diag(ones(4,1))*1E-10;
%             N = 0;
%             [K_lqr,~,~] = lqr(A,B,Q,R,N);
            
            % GET THE DESCRETE LQR CONTROLLER
            Q = diag([1 1 1 1 1 1 1 1 1 1 1 1])*1E1;
            R = diag(ones(4,1))*1E-10;
%             R = diag(ones(4,1))*1E-10;
            N = 0;
            [K_lqr,~,~] = dlqr(A,B,Q,R,N);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;psi;the;phi;v;u;w;p;q;r]



















