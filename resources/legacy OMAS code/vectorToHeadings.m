        % GET THE PITCH AND HEADING ERROR FROM THE INTERVAL INPUTS
        function [theta,lambda] = getControlInputs(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS
            % V - The current unit velocity vector (ie [1;0;0])
            % U - The desired heading vector       (ie midrad(0,[3,5,6])
            
            
            %% GET THE LINE OF SIGHT ANGLE             
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = diag([1 1 0])*V;
            Uh = diag([1 1 0])*U;                                          % Reject the vertical elements
            % DEFINE LOS ANGLE                        
            rotationAxis = cross(Vh,Uh);
            rotationSign = sign(mid(rotationAxis(3)));        
            % CALCULATE THE RELATIVE HEADING ANGLE
            lambda = rotationSign*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            if isnan(lambda)
                lambda = 0;
            end
            % GET THE ELEVATION ANGLE
            if isintval(U)
                theta = atan(U(3)/norm(Uh));
            else
                theta = atan2(U(3),norm(Uh));
            end
            if isnan(theta)
                theta = 0;
            end
        end