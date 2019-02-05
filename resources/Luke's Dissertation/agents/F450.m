% This class is designed as a child of the agent class "agent", used to 
% define the properties of the F450 quadrotor, as an agent for the purpose 
% of multi-vehicle control simulation.

% Author: James A. Douthwaite

% link: http://uk.mathworks.com/help/matlab/matlab_oop/class-constructor-methods.html

classdef F450 < agent
%%% F450 CHILD CLASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This class defines the F450 quadrotor specific properties for importing a
    % small quadrotor into the environment.
    
%% INITIALISE THE ADRONE SPECIFIC PARAMETERS
    properties
        % Phyiscal Parameters are initialised in the constructor
        ROTOR;
        % Performance Parameters
        battery_capacity    = 2500;     % Maximum battery capacity (mAh)
        battery_voltage     = 11.1;     % Maximum battery voltage (V)
    end
%%  CLASS METHODS
    methods 
        % CONSTRUCTION METHOD
        function obj = F450(namestr,mode)
            % Create AGENT with name tag, and incrementing ID number
             if nargin == 0
                namestr = ''; % Assign values to asset_args
             end
             %% INITIALISE THE AGENT SUPERCLASS
             obj@agent(namestr);
             %% INITIALISE THE BODY (PHYSICAL) PARAMETERS
             obj.BODY.armlength = (3.226E-01);			  % Measured directly from the CAD model (m)
             obj.BODY.width  = (obj.BODY.armlength)*2;	  % Assumed double the arm length
             obj.BODY.length = obj.BODY.width;
			 obj.BODY.mass = 1.2238;              % MTOW (kg)
             % Airframe Inertial Parameters
             switch mode
    			case '+'
        			% Inertial matrix (+) configuration
             		obj.BODY.I = 1.0E-06.*[2.7934E+04,1.5878E+01,3.9676E+00;
                               			   1.5878E+01,5.4873E+04,1.2336E+01;
                               			   3.9676E+00,1.2336E+01,2.8122E+04];
			    case 'x'
			        %Inertial matrix (x) configuration
			        obj.BODY.I = 1.0E-06.*[2.8032E+04,-1.5213E+01,9.3489E+01;
			                              -1.5213E+01,5.4973E+04,6.9167E+00;
			                               9.3489E+01,6.9167E+00,2.8023E+04]; 
			    otherwise
			        error('|| ERROR: Invalid F450 configuration (x : +).');
			end
			obj.BODY.I_inv = inv(obj.BODY.I);

			%% INITIALISE ROTOR PROPERTIES
			% Geometry
			obj.ROTOR.R = 1.2280E-01;	 % Rotor Radius (m)
			obj.ROTOR.A = 1.4000E-03;	 % Rotor Cross-section (m^2)
			obj.ROTOR.h = 5.8912E-02;	 % Rotor Elevation (m)
			obj.ROTOR.k = 1.4530E-05;	 % Thrust constant 1-5N
			obj.ROTOR.kt= 0.0529;		 % Linear region gradient
			obj.ROTOR.kh= 1.0000E-04;    % Constant: yaw moment- angular velocity

        end        
    end
%     enumeration
%         F450,drone,quadcopter,quadrotor
%     end
end