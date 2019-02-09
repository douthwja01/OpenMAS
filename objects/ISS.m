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
        % CONSTRUCTOR METHOD
        function obj = ISS(varargin)
            % This function is called to create the 'agent_example' class,
            % which then takes on the new parameters specified.
            
            % INPUT HANDLING
            if length(varargin) == 1 && iscell(varargin)                   % Catch nested cell array inputs
                varargin = varargin{:};
            end 
            % CALL THE SUPERCLASS CONSTRUCTOR
            obj@agent(varargin);                                           % Create the super class 'agent' 
            % GET THE DYNAMICS PROPERTIES
            [obj.DYNAMICS] = obj.getDynamicsProperties();
            
            % VIRTUAL DEFINITION
            obj.VIRTUAL.radius = sqrt((obj.length/2)^2 + (obj.width/2)^2 + (obj.height/2)^2);
            
            % CHECK FOR USER OVERRIDES
            [obj] = obj.configurationParser(obj,varargin);
            
            % SCALING
            [obj.GEOMETRY] = OMAS_graphics.scale(obj.GEOMETRY,[obj.length/2,obj.width/2,obj.height/2]); % To match real world dimensions
        end
        % AGENT MAIN CYCLE 
        function [obj] = main(obj,TIME,varargin)
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
            
            
            % GET THE TIMESTEP
            if isstruct(TIME)
                dt = TIME.dt;
            else
                error('Object TIME packet is invalid.');
            end
            % CHECK FOR NEW INFORMATION UPDATE
            observationSet = varargin{1}; % The detected objects                                
            if isempty(observationSet)
                % NO INFORMATION AVAILABLE
                passiveStateUpdate = obj.stateDynamics_accelerations(dt,[0;0;0],[0;0;0]); % Update the state vector
                obj = obj.updateGlobalProperties(dt,passiveStateUpdate);    % Update the objects global porperties
                return
            else
                % UPDATE THE AGENT WITH THE NEW INFORMATION
                [obj,~,~] = obj.getAgentUpdate(observationSet);
            end
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
            % INSERT ALGORITHM/DECISION MAKING PROCESS HERE
            orbitalAltitude = 6.371E+06 + 406E3; 
            orbitalTangentialVelocity = 7.67E3; % Orbital speed
            orbitalOmega = orbitalTangentialVelocity/orbitalAltitude;
            
            % AS STATE UPDATES
            velocity_k_plus = [0;orbitalTangentialVelocity;0];             % Moving sideways (assume facing the earth)
            omega_k_plus = [0;0;orbitalOmega];
            
            % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            
            % USE OUTPUT TO DEFINE NEW AGENT STATE
            [dXdt] = obj.dynamics_simple(obj.localState,velocity_k_plus,omega_k_plus);
            % SIMPLE INTEGRATION
            eulerState = obj.localState + dt*dXdt;
            
            % UPDATE THE 'agent_example' PROPERTIES WITH ITS NEW STATE
            [obj] = obj.updateGlobalProperties_ENU(dt,eulerState);
        end
    end
    
    methods (Static)
        % GET THE REPRESENTATIVE DYNAMICS PROPERTIES OF THE ISS
        function [DYNAMICS] = getDynamicsProperties() 
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