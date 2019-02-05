% This class is designed as a child of the agent class "agent", used to 
% define the properties of the AR quadrotor, as an agent for the purpose 
% of multi-vehicle control simulation.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef ARdrone < agent
%%% ARdrone CHILD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the ARdrone specific properties for importing a
    % small quadrotor into the environment.
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties
        % Phyiscal Parameters are initialised in the constructor
        % Performance Parameters
        battery_capacity    = 1000;     % Maximum battery capacity (mAh)
        battery_voltage     = 11.1;     % Maximum battery voltage (V)
        flight_time         = 720;      % Rated flight time (s)
        ultrasound_freq     = 40E3;     % Ultrasonic range-finder frequency (Hz)
        ultrasound_range    = 6;        % Ultrasonic range-fidner range (m)
        camera_freq         = 30;       % Camera capture frequency (fps)
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = ARdrone(namestr)
            % Create AGENT with name tag, and incrementing ID number
            if nargin == 0
               namestr = ''; % Assign values to asset_args
            end
            % INITIALISE THE AGENT SUPERCLASS
            obj@agent(namestr);

%%          THE BODY (PHYSICAL) PARAMETERS
%           obj.BODY.mass   = 0.366;               % MTOW (no hull)
%           obj.BODY.width  = 0.450;               % Physical width (m)
%           obj.BODY.length = 0.290;               % Physical length (m)
%           % WITH OUTDOOR HULL
%           obj.BODY.mass   = 0.400;               % MTOW (with outdoor shell)
%           obj.BODY.width  = 0.452;
%           obj.BODY.length = 0.452;
%           % WITH INDOOR HULL
            obj.BODY.mass   = 0.436;               % MTOW (with indoor shell)
            obj.BODY.width  = 0.515;
            obj.BODY.length = 0.515;
            % INERTIAL MATRIX (x) configuration (FICTIONAL I MATRIX)
            obj.BODY.I = 1.0E-06.*[2.8032E+04,0.0000E+00,0.0000E+00;
                                   0.0000E+00,2.8032E+04,0.0000E+00;
                                   0.0000E+00,0.0000E+00,2.8023E+04]; 
        end
        % STATE INITIALISATION HANDLED IN THE SUPER CLASS 'agent'
        
        % ARDRONE SPECIFIC FUNCTIONS
        % Battery check function
        
    end
%     enumeration
%         AR,drone,quadcopter,quadrotor
%     end
end