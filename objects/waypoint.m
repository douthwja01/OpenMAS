%% GENERIC WAYPOINT CLASS (waypoint.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic waypoint object, to be placed
% within the simulation. Upon achieving the waypoint, the waypoint object
% will be removed from the simulation. 

% Author: James A. Douthwaite 05/06/2017

classdef waypoint < objectDefinition
%%% WAYPOINT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    properties
        % Default ownership table
        ownership = struct(...
            'objectID',NaN,...
            'name','all',...
            'priority',0); 
        % Associated agent [objectID, name, priority]
    end
    %% ///////////////////////// MAIN METHODS /////////////////////////////
    methods 
        % Constructor
        function [this] = waypoint(varargin)
            % This function is to construct the waypoint object using the
            % object defintions held in the 'objectDefinition' base class.
            % We want to allow a waypoint to be called into existance and
            % made valid for any number of agent ID's (objectID's).
                        
            % Call the super class
            this = this@objectDefinition(varargin);
            
            % Allocate way-point defaults
            this = this.SetGLOBAL('type',OMAS_objectType.waypoint);
            this = this.SetGLOBAL('hitBoxType',OMAS_hitBoxType.spherical);
            this = this.SetGLOBAL('symbol','v');            
            this = this.SetGLOBAL('priorState',this.localState); 
                        
            % IF SIMPLE CONSTRUCTION.. RETURN NOW
            if length(varargin) < 1
                return
            end 
            
            % DEFAULT CONFIGURATION PARAMETERS
            priority = 0;
            agentSet = [];

            tmp = strncmpi(varargin,'owner',5);
            if any(tmp)
                agentSet = varargin{find(tmp) + 1};
            end 
            tmp = strncmpi(varargin,'priority',5);
            if any(tmp)
                priority = varargin{find(tmp) + 1};
            end                    

            % ALLOCATE ASSOCIATED AGENT PROPERTIES
            if ~isempty(agentSet)
                for allocatedAgent = 1:length(agentSet)
                    [this] = this.CreateAgentAssociation(agentSet{allocatedAgent},priority);
                end
            else
                % NO SPECIFIC AGENTS, ASSIGN PRIORITY TO THE GENERAL CASE
                this.ownership.priority = priority;
            end
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end      
    end
    %% /////////////////// WAYPOINT OWNERSHIP FUNCTIONS /////////////////////
    methods
        % CREATE ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [this] = CreateAgentAssociation(this,agentObject,priority)
            % This method is used to add an association between the
            % waypoint and a given agent object.
            
            % INPUT HANDLING
            if ~exist('priority','var')
               priority = this.ownership.priority; 
            end
            
            % Create a blank ownership table
            newOwnership = this.CreateOwnershipTable();
            newOwnership.objectID = agentObject.objectID;   % Assign a objectID
            newOwnership.name     = agentObject.name;       % Assign a name
            newOwnership.priority = priority;               % Assign a priority
            
            % OVERRIDE 'all' CASE                        
            if isnan(this.ownership(1).objectID)
                % If this is the first allocated object, redefine                            
                this.ownership(1) = newOwnership;            
                this.SetGLOBAL('colour',agentObject.GetGLOBAL('colour'));
                return
            end
            % Check if ID is already present on the table
            ownershipIDs = [this.ownership.objectID];
            IDindex = ownershipIDs(ownershipIDs == agentObject.objectID);
            if ~isempty(IDindex)
                warning('Unable to add agent to waypoint, agent already associated.');
                return
            else
                % Otherwise append the new ownership entry
                this.ownership = vertcat(this.ownership,newOwnership);       
            end
        end
        % REMOVE ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [this] = RemoveAgentAssociation(this,agentObject)
            % This function removes the waypoints association with a
            % specific objectID
            % INPUT HANDLING
            if nargin < 1
                error('Provide either an object ID or class object');
            elseif isnumeric(agentObject)
                agentID = agentObject;
            else
                agentID = agentObject.objectID;
            end
            
            % GET THE AGENTS ASSOCIATION
            [~,IDindex] = this.GetAgentAssociation(agentID);
            
            if isempty(IDindex)
                warning('Cannot remove; agent not associated with waypoint.');
                return
            end
            % REMOVE THE AGENT FROM THE WAYPOINT OWNERSHIP
            this.ownership(IDindex) = [];
        end   
        % FETCH ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [isAssociated] = IDAssociationCheck(this,agentID)
            % Simple binary check if is waypoint is associated with the
            % provided object ID number.           
            
            % Input sanity check
            assert(isnumeric(agentID),'Agent ID must be a valid objectID');
            % CHECK THE OWNERSHIP INDEX FOR THIS ID
            isAssociated = logical(any(this.ownership.objectID == agentID));
        end
        % FETCH ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [ownershipItem,IDindex] = GetAgentAssociation(this,agentObject)
            % This method is used to remove an association between the
            % waypoint and a given agent object.
            % INPUTS:
            % agentObject - Either the agent object or objectID.
            % OUTPUT:
            % this         - The updated waypoint object
            
            % INPUT HANDLING
            if nargin < 1
                error('Provide either an object ID or class object');
            elseif isnumeric(agentObject)
                agentID = agentObject;
            else
                agentID = agentObject.objectID;
            end
            % CREATE VECTOR OF ASSOCIATED IDs
            ownershipIDs = [this.ownership.objectID];
            % IF ASSOCIATION IS NOT THERE, THROW WARNING
            IDindex = find(ownershipIDs == agentID);
            if isempty(IDindex)
                warning('Agent is not associated with waypoint.');
                ownershipItem = [];
                return
            end
            % RETURN THE ASSOCIATION
            ownershipItem = this.ownership(IDindex);
        end
    end
    %% Private waypoint methods
    methods (Static, Access = private)
        % Create a new ownership table
        function [ownership] = CreateOwnershipTable()
            % Initialise the ownership table
            ownership = struct(...
                'objectID','',...
                'name','',...
                'priority',0);
        end
    end
end %%% STATE VECTOR IS DEFINED AS [x;y;z;v;u;w;psi;the;phi;p;q;r] %%%%%%%% 

