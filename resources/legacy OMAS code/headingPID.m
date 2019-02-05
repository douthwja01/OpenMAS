            %% \\\\\\\\\\\\\\ CALCULATE CONTROL INPUTS \\\\\\\\\\\\\\\\\\\\
            desiredHeading = desiredVelocity/norm(desiredVelocity);
            [e_theta,e_lambda] = obj.getControlInputs([1;0;0],desiredHeading);
            e_control = [0;e_theta;-e_lambda]; % -ve in the NED control frame
            % BUILD THE CONTROL INPUTS
            Kp = 0.2; Kd = 10;    % 0.2 , 10
            a_angular = Kp*eye(3)*e_control + Kd*eye(3)*(e_control-obj.priorError);
            % REMEMBER PREVIOUS ERROR
            obj.priorError = e_control;