
% This function provides the interface between the host system and the
% desired directory.

classdef OMAS_system
    
    properties
        % Default system description
        system = OMAS_systemType.windows;
        delimiter = '\';
        root;
    end
  
    methods (Static)
        % Get the system type enumeration
        function [systemEnum,obj] = GetSystemType(obj)
            % Determine the OS version
            switch 1
                case ispc
                    obj.system = OMAS_systemType.windows;
                    obj.delimiter = '\'; % DOS delimiter
                case ismac
                    obj.system = OMAS_systemType.OSX;
                    obj.delimiter = '/'; % Unix delimiter
                case isunix
                    obj.system = OMAS_systemType.linux;
                    obj.delimiter = '/'; % Unix delimiter
                otherwise
                    error('System architecture not recognised.');
            end
            systemEnum = obj.system;
            obj.rootPath = getenv('SYSTEMROOT');
        end
        % Get the system information
        function [archstr] = GetSystemInfo()
            archstr = computer('arch');
        end
        % Get the matlab installation information
        function [release,version,date] = GetMatlabInfo()
            installStruct = ver('matlab');
            release = installStruct.Release;    % Matlab release
            version = installStruct.Version;    % The release version
            date    = installStruct.Date;       % Date of the release
        end
    end
    % /////////////////////////// UTILITIES ///////////////////////////////
    methods (Static)
        % Get a the program dependancies
        function [isSuccessful] = GetFileDependancies()
            % This function ensures that OpenMAS has access to the complete
            % set of file dependancies.
            
            % Get the path to the install directory
            repoString = mfilename('fullpath');
            ind = strfind(repoString,'environment');
            repoString = repoString(1:(ind-1));
            
            % Known file dependancies
            relativePaths = {...
                'environment',...
                'objects',...
                'scenarios'...
                'toolboxes',...
                'environment\assets',...
                'environment\events',...
                'environment\assets',...
                'environment\common',...
                'environment\figures'};
            
            % Parse known matlab paths
            pathCell = regexp(path, pathsep, 'split');                                 % The name of the paths
            
            % Move through the relative path list
            for i = 1:numel(relativePaths)
                % The dependency paths
                dependencyPath = OMAS_system.GetOSPathString([repoString,relativePaths{i}]); 
                % If it the dependancy is not on the system path
                if any(strcmpi(dependencyPath, pathCell))
                    continue
                end
                % Attempt to add the dependancy
                try
                    addpath(dependencyPath);
                catch depError
                    warning(depError.message);
                    isSuccessful = 0;
                    return
                end
            end
            % Indicate successful
            isSuccessful = 1;
        end
        % Get a string path in the correct OS
        function [pathString] = GetOSPathString(pathString)
            % Input sanity check
            assert(ischar(pathString),'Provided string is assumed to be a path.');
            % Get the system type
            [~,sys] = OMAS_system.GetSystemType();
            % Format the path in the OS version                   
            pathString = strrep(pathString,'/',sys.delimiter);
            pathString = strrep(pathString,'\',sys.delimiter);
        end 
    end
end