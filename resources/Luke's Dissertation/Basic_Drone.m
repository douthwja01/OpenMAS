%%

classdef Basic_Drone < agent
    %%% Basic Drone Child Class %%%
        % Basic requirements for drone to be use in the
        % collision avoidance algorithm
    
    properties
        Size = 1;          % Size of drone at largest part in meters      
    end
    
    methods
        function obj = Basic_Drone(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
               namestr = ''; % Assign values to asset_args
            end
            % INITIALISE THE AGENT SUPERCLASS
            obj@agent(namestr);
            
        end
    
    end
end
