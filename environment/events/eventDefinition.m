%% EVENT OBJECT (event.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class is designed to define a generic agent and import this variables 
% into the simulation space for the purpose of multi-vehicle control simulation.

% Author: James A. Douthwaite

classdef eventDefinition
%%% EVENT BASE CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class contains the basic properties of a generic event, neither
    properties
        % EVENT PARAMETERS
        eventID;           
        name;
        type;
        time;
        info;
    end
%%  CLASS METHODS
    methods 
%       CONSTRUCTION METHOD
       function obj = eventDefinition(time,name,summaryInfo)
            % This function is designed to generate an event with an 
            % automatically defined ID number and basic properties
            % INPUTS:
            % time
            % name
            % summaryInfo - Assigned description string
            % OUTPUTS:
            % obj     - The generated object
            
            % INPUT HANDLING
            % CONFIRM EVENT TIME
            if ~isnumeric(time)                    % Default time setting
                %             if ~exist('time','var') || ~isnumeric(time)                    % Default time setting
                obj.time = 0;
                %warning('An event time must be defined.');
                return
            else
                obj.time = time;
            end
            % CONFIRM NAME STRING
            if isnumeric(name) || isempty(name)    % Default name setting
                obj.name = 'EVENT';
            else
                obj.name = name;
            end
            % CONFIRM DEFAULT SUMMARY INFORMATION
            %             if ~exist('summaryInfo','var') || isempty(summaryInfo)
            if isempty(summaryInfo)
                obj.info = 'None.';                                        % Default information setting
            else
                obj.info = summaryInfo;
            end
            
            % ALLOCATE AUTOINCREMENTING ID NUMBER
            persistent eventcount;
            if isempty(eventcount)
                eventcount = 1;
            else
                eventcount = eventcount + 1;
            end
            obj.eventID = uint16(eventcount);   % Assign the number as an ID tag
            
            % ASSIGN THE DEFAULT TYP
            obj.type = eventType.event;         % Assign default enumeration of 'event'
       end
    end
    % PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        % PRINT THE OBJECT SUMMARY
        function displaySummary(obj)
            fprintf('[%s]\tevent created (ID:%d @%.2fs):\t%s\n',obj.name,obj.eventID,obj.time,obj.info);
        end
    end
end

