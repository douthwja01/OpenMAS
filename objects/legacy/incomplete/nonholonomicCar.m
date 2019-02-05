%% NON-HOLONOMIC CART AGENT MODEL (nonholonomicCar.m) %%%%%%%%%%%%%%%%%%%%%
% This car 

% Author: James A. Douthwaite

classdef nonholonomicCar < agent
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties    

    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD FOR 
        function obj = nonholonomicCar(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
               namestr = '';                              % Assign default naming scheme
            end
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(namestr);                           % Get super class 'agent'
        end
        %% AGENT MAIN CYCLE 
        function [obj] = processTimeCycle(obj,TIME,varargin)
            % This function is designed to house a generic agent process
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >dt      - The timestep
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % INPUT HANDLING
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            % CHECK FOR ENV UPDATE
%             environmentalData = varargin{1};                                     % The detected objects          
%             % CALCULATE THE CONTROL LOOP CYCLE
%             if 0 == rem(TIME.currentTime,(1/obj.sampleFrequency)) 
%                 % WAIT FOR LOOP CYCLE SCHEDULE
%                 obj = updateStateVector(obj,dt,[0;0;0],[0;0;0]);
%                 return
%             end
%             % IF NO OBJECTS IN RANGE YET
%             if isempty(environmentalData)
%                  % NO DETECTIONS, CONTINUE ON CURRENT COURSE
%                 obj = updateStateVector(obj,dt,[0;0;0],[0;0;0]);
%                 return
%             end
            dt = 0.1
            % 8. UPDATE THE STATE VECTOR
            obj = updateStateVector(obj,dt);
        end
        
        function obj = updateStateVector(obj,dt)
            
            obj.
        end
    end
    
    %%  PRIVATE METHODS (CLASS SPECIFIC TOOLS)
    methods (Access = private)
        % DEFINE THE 2D NON-HOLONOMIC MOVEMENT
        function obj = 2DNonHolonomicDynamics(obj,inputs)
            
        end
    end
    % STATIC MATHEMATIC FUNCTIONS
    methods (Static)
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]