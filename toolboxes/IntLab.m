%% IntLab Toolbox wrapper (IntLab.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef IntLab < ToolBox
   
    methods
        % Constructor
        function [obj] = IntLab(varargin)
            % Check Intlab is present in the environment
            try
                test = infsup(0,1);
            catch
                warning('Intlab tools not initialised. Bring up now...');
                
                repoPath = obj.GetRepoPath();
                % Attempt to load Int-lab in its relative directory
                obj.GetIntlab(repoPath);
            end 
        end
    end
    % Initialisation of the toolbox
    methods

    end
    %% ////////////////////// AUXILLARY METHODS /////////////////////////// 
    methods (Static, Access = private)
        % Load the intlab library as a wrapper
        function [successFlag] = GetIntlab(path)            
            
            % Add system paths

            addpath(OMAS_system.GetOSPathString([path,'toolboxes\Intlab_V7.1']));

            % Try to initialise the library
            try 
                run('startup.m'); 
            catch 
                warning('Failed to execute "startup.m" script.');
            end
            
            % Test the library is functioning
            try 
                test = infsup(0,1);
                clearvars test;
                successFlag = 1;
            catch intervalError
                warning("Unable to load the Interval library IntLab");
               % rethrow(intervalError);
                successFlag = 0;
            end
        end
    end 
end