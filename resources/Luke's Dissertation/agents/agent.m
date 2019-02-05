% This class is designed to define a generic agent and import this variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef agent
%%% AGENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic agent, neither
    % aerial or ground based.
    
    % Define generic agent properties
    properties
        agentID;
        name;
        BODY;           % Physical Body Parameters
        state;
        position;
        attitude;       % Initialise position & attitude
        % accessed via obj.mass etc.
    end
%%  CLASS METHODS
    methods 
%       CONSTRUCTION METHOD
        function obj = agent(namestr)
%           Create AGENT with name tag, and incrementing ID number
            persistent agentcount;
            if isempty(agentcount)
                agentcount = 1;
            else
                agentcount = agentcount + 1;
            end
            obj.agentID = agentcount;   % Assign the number as an ID tag
            
%           No input name string specified, use greek naming scheme
            if nargin == 0 || strcmp(namestr,'')
                defaultID = { 'alpha', 'beta','gamma', 'delta','epsilon'  ...
                            ,  'zeta',  'eta','theta',  'iota','kappa'    ...
                            ,'lambda',   'mu',   'nu',    'xi','omicron'  ...
                            ,    'pi',  'rho','sigma',   'tau','upsilon'  ...
                            ,   'phi',  'chi',  'psi', 'omega'};
%               Generate default name string
                deflength = length(defaultID);
                cycle = ceil(agentcount/deflength)-1;
                IDval = agentcount - (cycle*deflength);         %  start from 'Alpha' and '001'
                namestr = strcat(defaultID{IDval},num2str(cycle+1,'%03d'));

%               fprintf('IDval: %d\t',IDval);
%               fprintf('Cycle: %d\t',cycle);
              fprintf('agentcount: %d\t',agentcount);
              fprintf('name: %s\n',namestr);
            end
            obj.name = namestr;
            
        end
%       INITIALISE STATE METHOD [X;Y;Z;sig;the;phi]
        function obj = initialise(obj,state)
            if 1 == size(state)
                state = transpose(state);
            end
            % Generate initial position vectors
            obj.position = state(1:3,1);
            obj.attitude = state(4:6,1);
            obj.state = vertcat(obj.position,obj.attitude); 
        end
%       UPDATE AGENT POSITION VECTOR
        function obj = update_position(obj,pos)
            % Check length of position vector
            if length(pos) == 3
                obj.position = pos;
            elseif length(pos) > 3
                % Assume the position is the first three values
                obj.position = pos(1:3,1);
            else
                error('|| Position vector unreadable.');
            end
        end
%       UPDATE AGENT ATTITUDE VECTOR
        function obj = update_attitude(obj,atti)
            % Check length of position vector
            if length(atti) == 3
                obj.attitude = atti;
            elseif length(atti) > 3
                % Assume the angles are given following position
                obj.attitude = atti(4:6,1);
            else
                error('|| Attitude vector unreadable.');
            end
        end
%       UPDATE AGENT STATE VECTOR
        function obj = update_state(obj,state)
            if length(state) >= 6 
                obj.position = state(1:3,1);
                obj.attitude = state(4:6,1);
                obj.state = vertcat(obj.position,obj.attitude);
            else
                error('|| State vector unreadable.');
            end
        end
        % Remove vehicle    
        function delete(obj)
             delete(obj)   
        end
    end
end

