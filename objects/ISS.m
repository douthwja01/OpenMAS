%% THE INTERNATIONAL SPACE STATION (ISS.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If have any questions or help getting started, contact me directly at:
% jadouthwaite1@sheffield.ac.uk

% Author: James A. Douthwaite

classdef ISS < agent
    properties
        orbit = 6.371E+06 + 406E3;      % Orbit above Earth
        orbitalSpeed = 7.67E3;          % Orbital speed
        length = 72.8;                  % Length (m)
        width = 108.5;                  % Width (m)
        height = 20;                    % Height (m)
        scale = 1;
    end

    methods 
        % Constructor
        function [this] = ISS(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
            
            % Call the super class
            this@agent(varargin);                                           
            
            % Assign defaults
            this.DYNAMICS = this.CreateDYNAMICS();
            this.GEOMETRY = OMAS_graphics.scale(this.GEOMETRY,[this.length/2,this.width/2,this.height/2]); % To match real world dimensions
            this.radius   = (sqrt((this.length/2)^2 + (this.width/2)^2 + (this.height/2)^2));
            
            % //////////////// Check for user overrides ///////////////////
            [this] = this.ApplyUserOverrides(varargin); % Recursive overrides
            % /////////////////////////////////////////////////////////////
        end 
        % Setup
        function [this] = setup(this,v,eta)
            this.localState = zeros(6,1);
            this.localState(4:6,1) = eta;
            % Define the prior state
            this.SetGLOBAL('priorState',this.localState);
        end
        % Main 
        function [this] = main(this,ENV,varargin)
            % This function is designed to contain everything your agent does
            % in a given simulation timestep. As an 'agent', a list of
            % detected obstacles is given if detected.
            % cycle that results in an acceleration vector in the global axis.
            % INPUTS:
            % varargin - Cell array of inputs
            % >TIME    - The current TIME structure (currentTime,dt, etc..)
            % >objects - The detectable objects cell array of structures
            % OUTPUTS:
            % obj      - The updated project
            
            % CHECK FOR NEW INFORMATION UPDATE
            [this,~,~] = this.GetAgentUpdate(ENV,varargin{1});
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
            % INSERT ALGORITHM/DECISION MAKING PROCESS HERE
            orbitalAltitude = 6.371E+06 + 406E3; 
            orbitalTangentialVelocity = 7.67E3; % Orbital speed
            orbitalOmega = orbitalTangentialVelocity/orbitalAltitude;
            
            % AS STATE UPDATES
            dt = ENV.dt;
            velocity_k_plus = [0;orbitalTangentialVelocity;0];             % Moving sideways (assume facing the earth)
            omega_k_plus    = [0;0;orbitalOmega];
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            X = this.localState(1:6);
            U = [velocity_k_plus;omega_k_plus];
            % USE OUTPUT TO DEFINE NEW AGENT STATE
            [dXdt] = this.SingleIntegratorDynamics(X,U);
            % SIMPLE INTEGRATION
            eulerState = this.localState + dt*dXdt;
            this.localState = eulerState;
            %this.localState(1:6)  = this.localState(1:6) + dt*dXdt;
            %this.localState(7:12) = dXdt;
            
            % UPDATE THE 'agent_example' PROPERTIES WITH ITS NEW STATE
            [this] = this.GlobalUpdate(dt,eulerState);
        end
    end
    methods (Static)
        % GET THE REPRESENTATIVE DYNAMICS PROPERTIES OF THE ISS
        function [DYNAMICS] = GetDynamicsProperties() 
            % Some general information:
%             orbitPerigee 403km
%             orbitApogee 406km
%             orbitInclination = deg2rad(51.64);      % Oribit inclination relative the earths lateral axis
%             orbitalSpeed = 7.67E3;                  % Orbital speed (km/s)              
%             orbitalPeriod = 92.49*60;               % Orbital period (s)

            % CONSTRUCT THE DYNAMICS CONTAINER
            DYNAMICS = struct('mass',419725);       % Total Mass (kg)            
        end
    end
end
% AGENT STATE VECTOR [x;y;z;v;u;w;psi;the;phi;p;q;r]