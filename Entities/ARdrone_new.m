
classdef ARdrone_new < quadcopter
    properties
        config = 'OUTDOOR'; % Indoor/outdoor dynamics
        % Performance Parameters
        battery_capacity  = 1000;	% Maximum battery capacity (mAh)
        battery_voltage   = 11.1;	% Maximum battery voltage (V)
        flight_time       = 720; 	% Rated flight time (s)
    end
    methods
        % Constructor
        function [this] = ARdrone_new(varargin)
            % Call the super class
            this@quadcopter(varargin);
            
            % Create the dynamics of a quadcopter
            this.DYNAMICS = this.CreateDYNAMICS();
            
            % //////////////// Check for user overrides ///////////////////
            this = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end
        %         % Setup
        %         function [this] = setup(this,localXYZVelocity,localXYZrotations)
        %         end
    end
    
    methods
        % Get the (generic) quadcopter dynamic properties
        function [DYNAMICS] = CreateDYNAMICS(this)
            
            % Import the default properties of a simple quadcopter
            DYNAMICS = this.CreateDYNAMICS_default();
            
            % ALTER THE DYNAMICS BASED ON THE CONFIGURATION
            switch upper(this.config)
                case 'NONE'
                    DYNAMICS.m = 0.366;
                case 'OUTDOOR'
                    DYNAMICS.m = 0.400; % Arm length equiv : 0.3196
                case 'INDOOR'
                    DYNAMICS.m = 0.436;
                otherwise
                    error('Configuration not recognised');
            end
            % Assign the ARdrone properties
            DYNAMICS.I = 1.0E-06.*[...
                2.8032E+04,0.0000E+00,0.0000E+00;
                0.0000E+00,2.8032E+04,0.0000E+00;
                0.0000E+00,0.0000E+00,2.8023E+04];
        end
    end
end