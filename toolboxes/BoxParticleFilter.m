%% Box Particle Filter (BoxParticleFilter.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class defines the behaviour of a simple Box Particle Filter (BPF)

classdef BoxParticleFilter < ParticleFilter & IntLab 
    % Properties of the box particle filter
    properties
        % Particle parameters -> defined in 'ParticleFilter'
        % n         - The number of particles
        % particles - The list of particle structures (x, z, w)
        
        % State/ Measurement Parameters
        % x0        - Initial state estimate
        % Qx;       % Covarience matrix (state update)
        z0;         % Initial measurement
        Vz;         % Initial measurement variance
        
        % BPF parameters
        dt = 0.1;           % Sampling rate

        
        q = 5e-2;           % State noise magnitude
        RANDOM_INOISE_ENABLE = true;
        n_factor = 1.0;     % The prediction noise factor   
        
        
%         % States
%         x_0;        % Initial state vector
%         x_k;
%         x_k_plus;   % [x y z dx dy dz]
%         % Measurements
%         h_k;        % [r dr theta phi]
%         % Interval measurements
%         z_k;        % [r_i r_s dr_i dr_s theta_i theta_s phi_i phi_s] 
    end
    
    methods
        % Constructor
        function [obj] = BoxParticleFilter(varargin)
            % The BPF for PF
            obj@ParticleFilter(varargin);
            
            % Pass the inputs to the input parser
            obj = obj.GetConfiguration(obj,varargin);
            
            % Initialise the BPF
            obj = obj.setup();
        end
        % Setup the BPF
        function [obj] = setup(obj)
            % This function computes the setup routine for the BPF.
            % x = [ x y z dx dy dz]
            
            obj.Vz = diag([2*pi/180;2*pi/180;5]);
            
            % State covariance matrix Qk-1
            Q1 = [obj.dt^3/3  obj.dt^2/2; 
                  obj.dt^2/2     obj.dt];
              
            % [TO-DO] Covarience aligned with states  
            obj.Qx = obj.q*kron(eye(3),Q1);   
            
            % Initialise the particle population
            obj.particles = obj.GetParticles_z(obj.z0);            
        end
        % Update - from new measurement
        function [obj] = update(obj,z_k)
            % The function cycles the box particle filter given a new
            % measurement at time k.
            
            % ////////////// UPDATE THE PARTICLE POPULATION ///////////////
            % Update the particle state
            [obj] = obj.ParticleStateUpdate();
            % Update the expected observations
            %[obj] = obj.ParticleMeasurementUpdate();
            % Update the weights
            [obj] = obj.ParticleWeightUpdate(z_k); % (with reference to new measurement)
            % /////////////////////////////////////////////////////////////
            
        end
    end
    
    %% ////////////////// BOX PARTICLE FILTER METHODS /////////////////////
    methods (Access = private)
        % Initialise the particles array (from measurement)
        function [particles] = GetParticles_z(obj,z0)
            % This function constructs a structure containing the set of
            % particles their states, and their weights in accordance to
            % the BPF.
            
            % Input sanity check
            assert(size(obj.Vz,2) == size(z0,1),...
                "The initial measurement varience is incompatable with the measurement vector.");
            
            % ////////////// DEFINE THE PARTICLE STRUCTURE ////////////////
            % The particle structure
            particles = struct(...
                'x_k',intval(),...      % The particles prior state
                'z_k',intval(),...      % The particle observation
                'w_k',1);               % The particle weight
            
            % Initialise the particle array
            for i = 1:obj.n
                % Generate random measurement
                zRand = z0 + sqrt(obj.Vz)*randn(size(z0));
                
                % Convert to initial particle measurement vector
                particles(i).z_k = midrad(zRand,diag(obj.Vz)/2); 

                % Estimate of initial state, from initial measurement
                xRandInt = obj.SphericalToCartesian(particles(i).z_k);
                
                % Assign the randomised particle states
                particles(i).x_k = [xRandInt;NaN(3,1)];    
                
                % Default weight
                particles(i).w_k = 1/obj.n;
            end
        end
        % Initialise the particles array (from state)
        function [particles] = GetParticles(obj,x0)
            % This function constructs a structure containing the set of
            % particles their states, and their weights.
            
            % Input sanity check
            assert(size(obj.Vx,2) == size(x0,1),"The initial state varience is incompatable with the state vector.");
            
            % The particle structure
            particles = struct(...
                'x_k',[],...      % The particle state
                'z_k',[],...      % The particle observation
                'w_k',[]);        % The particle weight
            
            % Initialise the particle array
            for i = 1:obj.n                
                % Generate random state
                xRand = x0 + sqrt(obj.Vx)*randn(size(x0));
                
                % Convert to initial particle measurement vector
                xRandInt = midrad(xRand,diag(obj.Vx)/2); 
                
                % Estimate of initial state from initial measurement
                particles(i).x_k = [xRandInt;NaN(3,1)];
                
                % Estimate of initial state, from initial measurement
                particles(i).z_k = obj.CartesianToSpherical(particles(i).x_k(1:3));
                
                % Default weight
                particles(i).w_k = 1/obj.n;
            end
        end
        % State-transition(prediction) model
        function [obj] = ParticleStateUpdate(obj)
            % This function determines the evolution of the state of a
            % given particle.
            
            for i = 1:obj.n 
                % The bounded noise              
                a = -obj.n_factor*sqrt(obj.Qx*eye(size(obj.Qx)));  
                b =  obj.n_factor*sqrt(obj.Qx*eye(size(obj.Qx)));   
                Inoise = infsup(a,b);
                % Generate the interval noise 
                Inoise = obj.RANDOM_INOISE_ENABLE*Inoise*randn(1);             % Choice of randomizing the interval noise or not

                % Define the state transistion matrix
                F = [1  0   0   obj.dt  0       0;
                     0  1   0   0       obj.dt  0; 
                     0  0   1   0       0       obj.dt; 
                     0  0   0   1       0       0;
                     0  0   0   0       1       0;
                     0  0   0   0       0       1];

                % Integrate the particle states as a prediction    
                obj.particles(i).x_k = F*obj.particles(i).x_k + Inoise;
            end
        end
        
%         % Observation model
%         function [obj] = ParticleMeasurementUpdate(obj)
%             % This function updates the particle measurements
%             for i = 1:obj.n
%                 obj.particles(i).z_k = obj.ObservationModel(obj.particles(i).x_k);
%             end
%         end
        
        
        % Particle weight update function
        function [obj] = ParticleWeightUpdate(obj,z_k)
            % This function recomputes the particles weightings based on
            % the BPF measurement contraction process.
            
            % Save the uncontracted 
            prev_particles = obj.particles;
            
            % Contraction loop
            for i = 1:8 
                % Use the measurement to contract the predictions
                obj = obj.MeasurementContraction(z_k);
            end
            
            % Assign new weights to the particles
            for i = 1:obj.n
                % If the state is nan, assign minimal weight
                if sum(isnan(obj.particles(i).x_k)) == numel(obj.particles(i).x_k)
                    obj.particles(i).w_k = 0;
                    continue
                end
                % Adjust the particle weighting
                obj.particles(i).w_k = obj.particles(i).w_k*prod(diam(obj.particles(i).x_k))/prod(diam(prev_particles(i).x_k));
            end
            % Re-normalise the particle weights
            obj.NormalizeParticleWeights();
        end
    end
    
    % BPF 
    methods
        % Measurement contraction process
        function [obj] = MeasurementContraction(obj,z_k)
            % This function completes the measurement contraction process
            % based on the current state predictions.
            
            
        end
    end
    
    % Supporting math functions
    methods (Static)
        % Divide interval 
        function [p_out] = subdivide_box(p,d,n)
            inf_d = inf(p(d));
            sup_d = sup(p(d));
            ints_d = linspace(inf_d,sup_d,n+1);
            p_out(d,1:n) = infsup(ints_d(1:end-1),ints_d(2:end));

            K=[1 2 3 4];
            K_minus_d = find(K-d); %for instance if k=2, K_minus_d=[1 3 4];
            p_out(K_minus_d,1:n)=repmat(p(K_minus_d),1,n);
        end
        % Divide interval (with zero handling)
        function [z] = divide_patch(y,x)
            % this function is a patch added to intlab for division operation.
            %in intlab, a division x/y and where y contains 0 can lead to NAN values
            % for instance
            %                 xx= infsup(0,2), yy= infsup(0,4), zz=xx/yy
            %                 intval xx =
            %                 [   0.00000000000000,   2.00000000000000]
            %                 intval yy =
            %                 [   0.00000000000000,   4.00000000000000]
            %                 intval zz =
            %                                 NaN
            
            % this is a correct spirit for intlab but this leads to incorrect values
            % for contraction algorithms.
            
            z = y./x;
            ii = find(inf(x).*sup(x)<=0);
            infinity = 999999999999;
            for i = ii
                if and(inf(x(i))<0,sup(x(i))>0)
                    % x(i) contains 0'
                    z(i)=infsup(-infinity, infinity);
                elseif and(inf(x(i))==0,sup(x(i))~=0)
                    dummy=infsup(1/sup(x(i)), infinity);
                    z(i)=y(i)*dummy;
                elseif and(inf(x(i))~=0,sup(x(i))==0)
                    dummy=infsup(-infinity, inf(x(i)));
                    z(i)=y(i)*dummy;
                end
            end
        end   
        % Convert Cartesian state to spherical state
        function [z_k] = CartesianToSpherical(xyz_k)

            r = norm(xyz_k);                                            % The range
            phi = atan2(xyz_k(2),xyz_k(1));                               % The angle made in the azimuth (bearing)
            theta = atan2(xyz_k(3),sqrt(xyz_k(1).^2 + xyz_k(2).^2));        % The elevation (vertical bearing)
            
            % For clarity
            z_k(1) = phi;
            z_k(2) = theta;
            z_k(3) = r;
            
            % Ensure column vector        
            if size(xyz_k,2) ~= size(z_k,2)
                z_k = z_k';
            end
        end
        % Convert spherical measurement to Cartesian state
        function [xyz_k] = SphericalToCartesian(z_k)
            % This function calculates the state estimate from the
            % measurement update.
            
            % For clarity
            psi_k   = z_k(1);
            theta_k = z_k(2);
            r_k     = z_k(3);
            
            % Get the cartesian state from circular measurement
            xyz_k = [cos(psi_k)*cos(theta_k);...
                     sin(psi_k)*cos(theta_k);...
                     sin(theta_k)]*r_k;
                    
            % Ensure column vector        
            if size(xyz_k,2) ~= size(z_k,2)
                xyz_k = xyz_k';
            end
        end
    end
end