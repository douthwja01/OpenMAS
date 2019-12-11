%% FORMATION CONTROLLED QUADCOPTER (quadcopter_formation.m) %%%%%%%%%%%%%%%
% This class is designed to implement the formation control elements of 
% 'agent_formation' with the dynamics defined by 'quadcopter'.

% Author: James A. Douthwaite

classdef quadcopter_formation < quadcopter & agent_formation
    properties
        % globalPosition - The position of the object in the global axes
        % globalVelocity - The velocity of the object in the global axes
        % quaternion     - The quaternion representing the earth to rotated body
        
        % DYNAMICS - All the models parameters are held in the DYNAMICS
        %            container field.
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = quadcopter_formation(varargin)
            
            % Call the super class
            this@quadcopter(varargin);                                      % Create the super class 'agent'

            % IMPORT THE AGENT DYNAMICS & CONTROL PARAMETERS
            this.DYNAMICS = this.CreateDYNAMICS();
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end    
        % Main
        function [this] = main(this,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
                        
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [this,~,agentSet,~] = this.GetAgentUpdate(ENV,varargin{1});                   % IDEAL INFORMATION UPDATE 
            
            if ~isempty(agentSet)
                [desiredVelocity,~] = this.formationControl_distance(agentSet);
            else
                desiredVelocity = [1;0;0];
            end
            
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [this] = this.controller(ENV,desiredVelocity);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            this.DATA.inputNames = {'$\dot{x}$ (m/s)','$\dot{y}$ (m/s)','$\dot{z}$ (m/s)',...
                                   '$\dot{\phi}$ (rad/s)','$\dot{\theta}$ (rad/s)','$\dot{\psi}$ (rad/s)'};
            this.DATA.inputs(1:numel(this.DATA.inputNames),ENV.currentStep) = this.localState(7:end);          % Record the control inputs 
        end
    end    
end