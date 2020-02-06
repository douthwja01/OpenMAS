classdef agent_3DRVO2 < agent_VO
    %% 3D ORCA properties
    properties
        
    end
    %% Main methods
    methods
        % Constructor
        function [this] = agent_3DRVO2(varargin)
            % Instantiate the parent agent
            this@agent_VO(varagin);
            
            
            [this] = ApplyUserOverrides(this,varargin);
        end
        % Main
        function [this] = main(this,varargin)
            
            
            
            
            % Pass the velocity to the controller
            this.controller(v_desired);
        end
    end
    %% Auxillary methods
    methods
        % The high-level avodiance rotation
        function [v_set] = GetAvoidanceTrajectory(this,v_pref)
            
            
        end
    end
end