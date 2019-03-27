%% EXAMPLE AGENT CLASS (agent_example.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This agent is intended to describe how to create a new agent class for
% simulation under the Open Multi-Agent Simulation (OpenMAS) framework. It 
% contains examples of the necessary properties and methods the class must
% have in order to interact with the environment.

% GENERAL NOTICE FOR AGENT DEFINITIONS
% - All objects(including agents) must be in '/objects' for dependancies.
% - The 'example_setup' demonstrates how the dependancies are setup.
% - All objects must have:
%   + A 'constructor' method (to create it).
%   + A 'main' method (to evaluate updates).
% - During the 'main' method, one cycle of the agents loop is computed,
%   updating the objects local state, and global properties. 

% To call the agent class, your work directory must be within '/objects',
% otherwise it must be added using 'addpath('objects')' to include the
% complete inheritance tree. 

% The class is called using the line:
% testAgent = agent_example();

% If have any questions or help getting started, please contact me directly 
% at: jadouthwaite1@sheffield.ac.uk

% Author: James A. Douthwaite

classdef agent_example < agent
    % SPECIFY YOUR AGENT-SPECIFIC PROPERTIES HERE
    % These are in addition to those inherited from 'agent'.
    properties
        % EXAMPLE PROPERTIES (Unique to this agent)
        % noOfWings = 2;
    end
    %% //////////////////// MAIN (required) METHODS ///////////////////////
    methods 
        % CONSTRUCTOR
        function obj = agent_example(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Create the super class 'agent'            
            % Set some of the parameters
            obj = obj.SetRadius(0.5);
            obj = obj.SetDetectionRadius(50);                              % Update the range attribute to the SIM VIRTUAL property            
            % Parse any overrides
            [obj] = obj.configurationParser(obj,varargin);
        end
        % AGENT MAIN CYCLE 
        function [obj] = main(obj,ENV,varargin)
            % This function is designed to contain everything your agent does
            % in a given simulation timestep. As an 'agent', a list of
            % detected entities is given if detected.
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV      - The simulation environmental update (currentTime,dt, etc..)
            % varargin - Cell array of optional inputs
            
            % OUTPUTS:
            % obj      - The updated project
            
            
            % GET THE TIMESTEP
            if isstruct(ENV)
                dt = ENV.dt;
            else
                error('OpenMAS environment structure is invalid.');
            end

            % PARSE THE ENTITY PROPERTIES FROM THE ENVIRONMENT 
            [obj,obstacleSet,agentSet] = obj.GetAgentUpdate(dt,varargin{1});       % IDEAL INFORMATION UPDATE       
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
            %        INSERT ALGORITHM/DECISION MAKING PROCESS HERE
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % DO NOTHING
            linearRates  = [0.5;0;0]; % (m/s)
            angularRates = [0;0;0.1]; % (rad/s)
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            %    DECIDE HOW THE AGENTS/ENTITIES LOCAL STATE IS UPDATED
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % How each agent updates its state can be defined locally, as
            % a dedicated function.
            % CALL THE SEPERATE UPDATE FUNCTION
            [state_k_plus] = obj.updateLocalState(ENV,...
                                       obj.localState,...
                                          linearRates,...
                                         angularRates);
                                     
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            %  USE THE UPDATED LOCAL STATE TO UPDATE THE GLOBAL PROPERTIES
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % How the objects global state is updated is dependant on how
            % the 'state_k_plus' is defined.
            obj = obj.updateGlobalProperties(dt,state_k_plus);
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % ASSIGNED STATE UPDATE FUNCTION (USING ODE45)
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