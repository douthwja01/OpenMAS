%% RECIPROCAL VO/FORMATION CONTROL AGENT (agent_formation_RVO.m) %%%%%%%%%%
% This agent is a hybrid of the formation control element agent class and
% the Reciprocal Velocity Obstacle (VO) class. 

% Author: James A. Douthwaite

classdef agent_formation_RVO < agent_formation_VO & agent_RVO
    %% INITIALISE THE AGENT SPECIFIC PARAMETERS
    properties
        
    end
    %%  CLASS METHODS
    methods
        % CONSTRUCTION METHOD
        function obj = agent_formation_RVO(varargin)
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@agent_formation_VO(varargin);
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
        end
        % NOTES:
        % This class simply overrides the 'process_timeCycle' function of
        % the superclass 'agent_formation_VO' such that it considers the
        % 'avoidanceCorrection' of the second super class 'agent_RVO'.
    end
end
% AGENT STATE VECTOR [x;y;z;phi;theta;psi;xdot;ydot;zdot;phidot;thetadot;psidot]