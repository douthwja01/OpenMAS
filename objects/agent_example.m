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
    
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = agent_example(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
                        
            % CALL THE SUPERCLASS CONSTRUCTOR
            this@agent(varargin);             % Create the super class 'agent'            
            
            % Set some of the parameters
            this.radius = 0.5;
            this.detectionRadius = 50;        % Update the range attribute to the SIM VIRTUAL property            
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        % Setup
        function [this] = setup(this,v,eta)
            % Define the initial state
            this.localState = zeros(6,1);
            this.localState(4:6,1) = eta;
        end
        % Main 
        function [this] = main(this,ENV,varargin)
            % This function is designed to contain everything your agent does
            % in a given simulation timestep. As an 'agent', a list of
            % detected entities is given if detected.
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % ENV      - The simulation environmental update (currentTime,dt, etc..)
            % varargin - Cell array of optional inputs
            
            % OUTPUTS:
            % obj      - The updated project
            
            % PARSE THE ENTITY PROPERTIES FROM THE ENVIRONMENT 
            [this,obstacleSet,agentSet] = this.GetAgentUpdate(ENV,varargin{1});       % IDEAL INFORMATION UPDATE       
            
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
            [state_k_plus] = this.UpdateLocalState(...
                ENV,...
                this.localState,...
                linearRates,...
                angularRates);
                                     
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            %  USE THE UPDATED LOCAL STATE TO UPDATE THE GLOBAL PROPERTIES
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            % How the objects global state is updated is dependant on how
            % the 'state_k_plus' is defined.
            this = this.GlobalUpdate(ENV.dt,state_k_plus);
            
            % ////////////// RECORD THE AGENT-SIDE DATA ///////////////////
            this.DATA.inputNames = {'$v_x$ (m/s)','$v_y$ (m/s)','$\dot{\psi}$ (rad/s)'};
            this.DATA.inputs(1:length(this.DATA.inputNames),ENV.currentStep) = this.localState(4:6);         % Record the control inputs
        end
    end
    %% /////////////////////// AUXILLARY METHODS //////////////////////////
    methods
        % ASSIGNED STATE UPDATE FUNCTION (USING ODE45)
        function [dXdt] = UpdateLocalState(this,TIME,X0,velocity,omega)
            % This function computes the state update for the current agent
            % using the ode45 function.
            
            U = [velocity;omega];
            dXdt = X0; % The default returned state
            
            % DETERMINE THE INTEGRATION PERIOD
            if TIME.currentTime == TIME.timeVector(end)
                return
            else
                [~,Xset] = ode45(@(t,X) this.SingleIntegratorDynamics(X,U),...
                    [0 TIME.dt],X0,...
                    odeset('RelTol',1e-2,'AbsTol',TIME.dt*1E-2));
                dXdt = Xset(end,:)'; % Pass the state at the end of the sample period
            end
        end
    end
end