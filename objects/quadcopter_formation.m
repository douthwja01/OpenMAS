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
    %  CLASS METHODS
    methods
        % CONSTRUCTOR METHOD
        function obj = quadcopter_formation(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@quadcopter(varargin);                                      % Create the super class 'agent'
            % INPUT HANDLING (Clean up nested loops)
            [varargin] = obj.inputHandler(varargin);
            % IMPORT THE AGENT DYNAMICS & CONTROL PARAMETERS
            [obj.DYNAMICS] = obj.importModelProperties();
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % //////////////////// AGENT MAIN CYCLE ///////////////////////////
        function [obj] = main(obj,ENV,varargin)
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % GET THE TIMESTEP
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('Object TIME packet is invalid.');
            end
            
            % //////////// CHECK FOR NEW INFORMATION UPDATE ///////////////
            % UPDATE THE AGENT WITH THE NEW ENVIRONMENTAL INFORMATION
            [obj,~,agentSet,~] = obj.GetAgentUpdate(ENV,varargin{1});                   % IDEAL INFORMATION UPDATE 
            
            if ~isempty(agentSet)
                [desiredVelocity,L] = obj.formationControl_distance(agentSet);
            else
                desiredVelocity = [1;0;0];
            end
            
            % PASS THE DESIRED VELOCITY TO THE DEFAULT CONTROLLER
            [obj] = obj.controller(dt,desiredVelocity);
            
            % /////////////// RECORD THE AGENT-SIDE DATA //////////////////
            obj.DATA.inputNames = {'u (m/s)','v (m/s)','w (m/s)','p (rad/s)','q (rad/s)','r (rad/s)'};
            obj.DATA.inputs(1:length(obj.DATA.inputNames),ENV.currentStep) = newState(7:12);         % Record the control inputs
        end
    end    
end