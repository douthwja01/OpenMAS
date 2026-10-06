%% 3D AGENT (agent_3D.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial agents. Planar agents live on agent_2D. Shared sensing, memory
% and waypoint handling stay on agent.

% Author: James A. Douthwaite

classdef agent_3D < agent
    methods
        function [this] = agent_3D(varargin)
            this@agent(varargin);

            this.localState = zeros(12,1);
            this.SetGLOBAL('is3D',true);
            this.DYNAMICS = this.CreateDYNAMICS();
            this = this.ApplyUserOverrides(varargin);
        end

        % Track a desired 3D velocity with the simple kinematic model.
        function [this] = Controller(this,dt,desiredVelocity)
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector is not a 3D numeric vector');
            assert(~any(isnan(desiredVelocity)),'The desired velocity vector contains NaNs.');

            [heading,speed] = this.nullVelocityCheck(desiredVelocity);
            [dPsi,dTheta] = this.GetVectorHeadingAngles([1;0;0],heading);
            omega = [0;dTheta;-dPsi]/dt;
            [omega_actual,speed_actual] = this.ApplyKinematicContraints(dt,omega,speed);

            if this.IsIdle()
                omega_actual = zeros(3,1);
                speed_actual = 0;
                this.v_nominal = 0;
            end

            [dX] = this.SimpleDynamics(this.localState(1:6),[speed_actual;0;0],omega_actual);
            this.localState(1:6)  = this.localState(1:6) + dt*dX;
            this.localState(7:12) = dX;
            this = this.GlobalUpdate_3DVelocities(dt,this.localState);
        end

        % Speed and heading PID controller.
        function [this] = Controller_PID(this,dt,desiredVelocity)
            assert(isnumeric(dt) && numel(dt) == 1,'The time step must be a numeric scalar.');
            assert(isnumeric(desiredVelocity) && numel(desiredVelocity) == 3,'Requested velocity vector must be a 3D local vector.');
            assert(~any(isnan(desiredVelocity)),'The requested velocity contains NaNs.');

            [unitDirection,desiredSpeed] = this.nullVelocityCheck(desiredVelocity);
            if abs(desiredSpeed) > this.v_max
                desiredSpeed = sign(desiredSpeed)*this.v_max;
            end

            [dPsi,dTheta] = this.GetVectorHeadingAngles([1;0;0],unitDirection);
            dHeading = [0;dTheta;-dPsi];
            e_speed = desiredSpeed - norm(this.localState(7:9));
            controlError = [e_speed;dHeading];

            Kp_linear = 0.8;
            Kd_linear = 0;
            Kp_angular = 1;
            Kd_angular = 0;
            control_fb = diag([Kp_linear Kp_angular Kp_angular Kp_angular])*controlError + ...
                diag([Kd_linear Kd_angular Kd_angular Kd_angular])*(controlError - this.priorError);
            this.priorError = controlError;

            speedFeedback = control_fb(1);
            omega = control_fb(2:4)/dt;
            [dX] = this.SimpleDynamics(this.localState(1:6),[speedFeedback;0;0],omega);
            this.localState(1:6)  = this.localState(1:6) + dt*dX;
            this.localState(7:12) = dX;

            if this.IsIdle()
                this.localState(7:12) = zeros(6,1);
            end
            this = this.GlobalUpdate_3DVelocities(dt,this.localState);
        end
    end
end
