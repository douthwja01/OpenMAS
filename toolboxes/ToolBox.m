%% Generic tools for integrating toolboxes with OpenMAS (toolbox.m) %%%%%%%

classdef ToolBox
    % Generic properties for toolbox definitions
    properties (Access = private)
        repositoryDir; 
    end
    
    methods
        % Constructor
        function [obj] = ToolBox(varargin)
            % Get the critical paths
            obj.repositoryDir = obj.GetRepoPath();            
            % Add the 'environment/OMAS' functions
            addpath([obj.repositoryDir,'environment']); 
            % Add the path to the 'common' tools
            addpath(OMAS_system.GetOSPathString([obj.repositoryDir,'environment\common']));   
            % Parse inputs against the object now the paths are added
            obj = obj.GetConfiguration(obj,varargin);
        end
    end
    %% ////////////////////// AUXILLARY METHODS ///////////////////////////
    methods (Static)
        % Parse the user inputs against the toolbox definition
        function [config] = GetConfiguration(defaultConfig,inputParameters)
            % This function is designed to parse a generic set of user
            % inputs and allow them to be compared to a default input
            % structure. This should be called in the class constructor
            
            % Input sanity check #1
            if nargin < 2 || isempty(inputParameters)
                config = defaultConfig;
                return
            end
            
            % Call the all parameter overrider
            try
                config = GetParameterOverrides(defaultConfig,inputParameters);
            catch ex
                warning('Has the common tools directory been pathed correctly?');
                rethrow(ex);
            end
        end
        % Get the relative paths
        function [masterPath] = GetRepoPath()
            % Get the system paths
            masterPath = mfilename('fullpath');
            ind = strfind(masterPath,'toolboxes');
            masterPath = masterPath(1:(ind-1));
        end
    end   
end