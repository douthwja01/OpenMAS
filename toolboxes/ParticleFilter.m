%% Particle Filter (ParticleFilter.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ParticleFilter < ToolBox
    
    properties
        % Particle parameters
        n = 4;      % Number of particles
        particles;  % Particle container
        
        % Parameters 
        x0;         % Initial state
        Vx;         % Initial state varience
        Qx;         % System noise covariance (in the state update)
        Qz;         % Measurement noise covariance (in the measurement)
        x_est;      % The state estimate
        
        % Other parameters
        k = 1;      % Measurement counter 
    end
    
    methods
        % Constructor
        function [obj] = ParticleFilter(varargin)
            % Get the parent class 
            obj@ToolBox(varargin)
            
            % Assign some default parameters
            obj.Vx = 2; % Initial state varience
            obj.Qx = 1; % System noise covariance (in the state update)
            obj.Qz = 1; % Measurement noise covariance (in the measurement)
            
            % Get the input configuration
            obj = obj.GetConfiguration(obj,varargin); 
            
            % Call the setup function
            obj = obj.setup();
        end        
        % Setup the PF
        function [obj] = setup(obj)
            % This function sets up the filter 
            
            % Input sanity check
            assert(numel(obj.x0) > 0 && isnumeric(obj.x0),'Please provide an initial state.');
            
            % Initialise the particles from the first state
            obj.particles = obj.GetParticles(obj.x0);
        end
        % Update - from new measurement
        function [obj] = update(obj,z_k)
            % This function cycles the particle filter given a new
            % measurement at time k. 
                        
            % ////////////// UPDATE THE PARTICLE POPULATION ///////////////
            % Update the particle state
            [obj] = obj.ParticleStateUpdate();
            % Update the expected observations
            [obj] = obj.ParticleMeasurementUpdate();
            % Update the weights
            [obj] = obj.ParticleWeightUpdate(z_k); % (with reference to new measurement)
            % /////////////////////////////////////////////////////////////
            
            % Randomly resample the particle states
            obj = obj.RandomResampling();
            % Estimate the state 
            obj.x_est(obj.k) = obj.StateEstimateUpdate();
            
            % Iterate the filter
            obj.k = obj.k + 1;
        end        
    end
    
    %% /////////////////// PARTICLE FILTER METHODS ////////////////////////
    methods (Access = private)
        % Initialise the particle array
        function [particles] = GetParticles(obj,x0)
            % This function constructs a structure containing the set of
            % particles their states, and their weights.
            
            % Input sanity check
            assert(size(obj.Vx,2) == size(x0,1),...
                "The initial state varience is incompatable with the state vector.");
            
            % The particle structure
            particles = struct(...
                'x_k',[],...      % The particle state
                'z_k',[],...      % The particle observation
                'w_k',[]);        % The particle weight
            
            % Initialise the particle array
            for i = 1:obj.n
                % Define the particle state
                particles(i).x_k = x0 + sqrt(obj.Vx)*randn(size(x0));
                % Default weight
                particles(i).w_k = 1;
            end
        end
        % State-transition model
        function [obj] = ParticleStateUpdate(obj)
            % The particle update is 
        	for i = 1:obj.n
                obj.particles(i).x_k = obj.StateTransitionModel(...
                    obj.k,...
                    obj.particles(i).x_k,...
                    obj.Qx);
            end
        end
        % Observation model
        function [obj] = ParticleMeasurementUpdate(obj)
            % This function updates the particle measurements
            for i = 1:obj.n
                obj.particles(i).z_k = obj.ObservationModel(obj.particles(i).x_k);
            end
        end
        % Particle weight update function
        function [obj] = ParticleWeightUpdate(obj,z_k)
            % This function updates the weighting of a given particle, in
            % relation to the new measurement z_k.
            
            % Assign gaussian weighting of the particles
            for i = 1:obj.n
                obj.particles(i).w_k = obj.GaussianWeighting(obj.Qz,z_k,obj.particles(i).z_k); 
            end
            % Renormalise the particle weightings
            [obj] = obj.NormalizeParticleWeights();
        end
        % Resampling
        function [obj] = Resample(obj)
            obj = obj.RandomResampling();
        end
        % State estimation model
        function [x_k_est] = StateEstimateUpdate(obj)
            % Simply use the mean of the particle states
            x_k_est = mean([obj.particles.x_k]); 
        end
    end
    
    % Utilities
    methods 
        % Randomly sample particles
        function [obj] = RandomResampling(obj)
            for i = 1:obj.n
                % Randomly sample the state
                randInd = find(rand <= cumsum([obj.particles.w_k]),1);
                % Modify the state of the particle to the new random state
                obj.particles(i).x_k = obj.particles(randInd).x_k;
            end
        end
        % Normalise the particle weights
        function [obj] = NormalizeParticleWeights(obj)
            % Normalise the new particle weights
            w_k_plus = [obj.particles.w_k]./sum([obj.particles.w_k]);
            for i = 1:obj.n
                obj.particles(i).w_k = w_k_plus(i);
            end    
        end
    end
    
    % Supporting math functions
    methods (Static)
        % Observation model (example)
        function [z_k_plus] = ObservationModel(x_k_plus)
            % This function defines the expected observation based on the
            % updates state.
            z_k_plus = x_k_plus^2/20;  
        end
        % State-transition model (example)
        function [x_k_plus] = StateTransitionModel(t,x_k,x_n)
            % This function provides an example state evolution for the
            % particles, based on the time, prior state and state varience.
            x_k_plus = 0.5*x_k + 25*x_k/(1 + x_k^2) + 8*cos(1.2*(t-1)) + sqrt(x_n)*randn;
        end
        % Weights defined by a gaussian equation
        function [w] = GaussianWeighting(Qz,z,z_update)
            w = (1/sqrt(2*pi*Qz)) * exp(-(z - z_update)^2/(2*Qz));
        end
    end
end