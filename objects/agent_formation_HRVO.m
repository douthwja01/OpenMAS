%% HYBRID-RVO/FORMATION CONTROL AGENT (agent_formation_HRVO.m) %%%%%%%%%%%%
% This agent is a hybrid of the formation control element agent class and
% the Hybrid-Reciprocal Velocity Obstacle (VO) class. 

% Author: James A. Douthwaite

classdef agent_formation_HRVO < agent_formation_RVO & agent_HRVO
    %% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods
        % CONSTRUCTION METHOD
        function obj = agent_formation_HRVO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_formation_RVO(varargin);
            
            % //////////////// Check for user overrides ///////////////////
            [obj] = obj.ApplyUserOverrides(varargin); % Recursive overrides
            % ///////////////////////////////////////////////////////////// 
        end
        % NOTES:
        % This class simply overrides the 'process_timeCycle' function of
        % the superclass 'agent_formation_RVO' such that it considers the
        % 'avoidanceCorrection' of the second super class 'agent_HRVO'.
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]