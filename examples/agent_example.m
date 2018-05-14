%% EXAMPLE AGENT CLASS (agent_example.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This agent is intended to describe how to create a new agent class for
% simulation under the Multi-Agent System Simulation (MASS) framework. It 
% contains examples of the necessary properties and methods the class must
% have in order to interact with the environment.

% To call the agent class, your work directory must be within '/objects',
% otherwise it must be added using 'addpath('objects')' to include the
% whole inheritance tree. 

% The class is called using the line:
% testAgent = agent_example

% Its associated properties and methods are accessed by:
% limitOfVision = testAgent.sensorRange
% updated_testAgent = testAgent.processTimeCycle(inputs)

% If have any questions or help getting started, contact me directly at:
% jadouthwaite1@sheffield.ac.uk

% Author: James A. Douthwaite

classdef agent_example < agent
    % SPECIFY YOUR AGENT-SPECIFIC PROPERTIES HERE
    % These are in addition to those inherited from 'agent'.
    properties
        % TEST PROPERTY
        % noOfWings = 2;
    end

    methods 
        %% CONSTRUCTOR METHOD
        function obj = agent_example(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Create the super class 'agent'            
            % AGENT SPECIFIC PARAMETERS
            % VIRTUAL DEFINITION
            obj.VIRTUAL.detectionRange = 50;                               % Update the range attribute to the SIM VIRTUAL property            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        %% AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to contain everything your agent does
            % in a given simulation timestep. As an 'agent', a list of
            % detected obstacles is given if detected.
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >TIME    - The current TIME structure (currentTime,dt, etc..)
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            % CHECK FOR NEW INFORMATION UPDATE
            observationSet = varargin{1}; % The detected objects                                
            if isempty(observationSet)
                % NO INFORMATION AVAILABLE
                passiveStateUpdate = obj.stateDynamics_accelerations(dt,[0;0;0],[0;0;0]); % Update the state vector
                obj = obj.updateGlobalProperties(dt,passiveStateUpdate);    % Update the objects global porperties
                return
            else
                % UPDATE THE AGENT WITH THE NEW INFORMATION
                [obj,~,~] = obj.getAgentUpdate(observationSet);
            end
            
            %% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
            % INSERT ALGORITHM/DECISION MAKING PROCESS HERE
            
            % DO NOTHING
            linearRates = [0.5;0;0];
            headingRates = [0;0;0.1];
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            
            % USE OUTPUT TO DEFINE NEW AGENT STATE
            state_k_plus = stateDynamics_singleIntegrator(obj,dt,linearRates,headingRates);
            
            % UPDATE THE 'agent_example' PROPERTIES WITH ITS NEW STATE
            obj = obj.updateGlobalProperties(dt,state_k_plus);
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]