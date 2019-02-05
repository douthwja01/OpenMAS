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
        
        % WAYPOINT PROPERTIES
        ownership = struct('objectID',NaN,...
                           'name','all',...
                           'priority',0); % Associated agent [objectID, name, priority]
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = waypoint(varargin)
            % This function is to construct the waypoint object using the
            % object defintions held in the 'objectDefinition' base class.
            % We want to allow a waypoint to be called into existance and
            % made valid for any number of agent ID's (objectID's).
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj = obj@objectDefinition(varargin); % Call the super class
            
            % WAYPOINT VIRTUAL DEFINITION
            obj.VIRTUAL.type       = OMAS_objectType.waypoint;            
            obj.VIRTUAL.hitBoxType = OMAS_hitBoxType.spherical;
            obj.VIRTUAL.symbol     = 'v';
            
            % IMPORT THE GEOMETRY IF IT EXISTS
            [obj.GEOMETRY] = obj.getObjectGeometry(obj);
      
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
                    [obj] = obj.createAgentAssociation(agentSet{allocatedAgent},priority);
                end
            else
                % NO SPECIFIC AGENTS, ASSIGN PRIORITY TO THE GENERAL CASE
                obj.ownership.priority = priority;
            end
            
            % CHECK FOR USER OVERRIDES
            obj.VIRTUAL = obj.configurationParser(obj.VIRTUAL,varargin); 
            obj = obj.configurationParser(obj,varargin);
        end           
    end
    % /////////////////// WAYPOINT OWNERSHIP FUNCTIONS /////////////////////
    methods
        % CREATE ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [obj] = createAgentAssociation(obj,agentObject,priority)
            % This method is used to add an association between the
            % waypoint and a given agent object.
            
            % INPUT HANDLING
            if ~exist('priority','var')
               priority = obj.ownership.priority; 
            end
            
            % CREATE NEW OWNERSHIP STRUCTURE
            newOwnership = struct('objectID',agentObject.objectID,...
                                      'name',agentObject.name,...
                                  'priority',priority);
            % OVERRIDE 'all' CASE                        
            if isnan(obj.ownership(1).objectID)                            
                obj.ownership(1) = newOwnership;    % If this is the first allocated object, redefine
                obj.VIRTUAL.colour = agentObject.VIRTUAL.colour;
                return
            end
            % CREATE VECTOR OF ASSOCIATED IDs
            ownershipIDs = [obj.ownership.objectID];
            % IF ASSOCIATION IS ALREADY THERE, DO NOT ADD
            IDindex = ownershipIDs(ownershipIDs == agentObject.objectID);
            if ~isempty(IDindex)
                warning('Unable to add agent to waypoint, agent already associated,');
                return
            else 
                obj.ownership = vertcat(obj.ownership,newOwnership);       % Otherwise append the new ownership
            end
        end
        % FETCH ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [isAssociated] = isIDAssociated(obj,agentID)
            % Simple binary check if is waypoint is associated with the
            % provided object ID number.           
            
            assert(isnumeric(agentID),'Agent ID must be a valid objectID');
            
            % CHECK THE OWNERSHIP INDEX FOR THIS ID
            isAssociated = any(obj.ownership.objectID == agentID);
        end
        % FETCH ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [ownershipItem,IDindex] = getAgentAssociation(obj,agentObject)
            % This method is used to remove an association between the
            % waypoint and a given agent object.
            % INPUTS:
            % agentObject - Either the agent object or objectID.
            % OUTPUT:
            % obj         - The updated waypoint object
            
            % INPUT HANDLING
            if nargin < 1
                error('Provide either an object ID or class object');
            elseif isnumeric(agentObject)
                agentID = agentObject;
            else
                agentID = agentObject.objectID;
            end
            % CREATE VECTOR OF ASSOCIATED IDs
            ownershipIDs = [obj.ownership.objectID];
            % IF ASSOCIATION IS NOT THERE, THROW WARNING
            IDindex = find(ownershipIDs == agentID);
            if isempty(IDindex)
                warning('Agent is not associated with waypoint.');
                ownershipItem = [];
                return
            else
                % RETURN THE ASSOCIATION
                ownershipItem = obj.ownership(IDindex);
            end
        end
        % REMOVE ASSOCIATION BETWEEN AGENT AND WAYPOINT
        function [obj] = removeAgentAssociation(obj,agentObject)
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
            [~,IDindex] = obj.getAgentAssociation(agentID);
            
            if isempty(IDindex)
                warning('Cannot remove; agent not associated with waypoint.');
                return
            else
                % REMOVE THE AGENT FROM THE WAYPOINT OWNERSHIP
                obj.ownership(IDindex) = [];
            end
        end
    end
end %%% STATE VECTOR IS DEFINED AS [x;y;z;v;u;w;psi;the;phi;p;q;r] %%%%%%%% 

