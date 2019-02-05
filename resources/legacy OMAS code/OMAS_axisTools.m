%% OPENMAS 3D-AXIS MANIPULATION TOOLS (OMAS_geometry.m) %%%%%%%%%%%%%%%%%%%
% This class provides the axis manipulation utilities used by OpenMAS to
% evaluate object progression over time.

% Author: James A. Douthwaite 06/09/2017

classdef OMAS_geometry
    % This class is designed to handle all the axis transformations and
    % conversions as part of the OpenMAS simulator.
    
    %% GENERIC TOOLS
    methods (Static) 
        % PLANAR ANGLE BETWEEN TWO VECTORS
        function [planarAngle] = planarAngle(u,v)
            planarAngle = acos(dot(u,v)/(norm(u)*norm(v)));
        end
        % ROTATE A VECTOR ABOUT AN AXIS VECTOR (RODRIGUES)
        function [v_rotated] = rodriguesRotation(u,k,theta)
            % v - Vector to be rotated
            % k - Is the rotation axis
            % Theta - The angle the vector is to be rotated through
            
            [m,n] = size(u);
            if (m ~= 3 && n ~= 3)
                error('input vector is/are not three dimensional')
            end
            if (size(u) ~= size(k))
                error('rotation vector v and axis k have different dimensions')
            end
            
            k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
            No = numel(u)/3; % number of vectors in array
            v_rotated = u; % initialize rotated vector array
            if ( n == 3 )
                crosskv = u(1,:); % initialize cross product k and v with right dim.
                for i = 1:No
                    crosskv(1) = k(2)*u(i,3) - k(3)*u(i,2);
                    crosskv(2) = k(3)*u(i,1) - k(1)*u(i,3);
                    crosskv(3) = k(1)*u(i,2) - k(2)*u(i,1);
                    v_rotated(i,:) = cos(theta)*u(i,:) + (crosskv)*sin(theta)...
                        + k*(dot(k,u(i,:)))*(1 - cos(theta));
                end
            else % if m == 3 && n ~= 3
                crosskv = u(:,1); % initialize cross product k and v with right dim.
                for i = 1:No
                    crosskv(1) = k(2)*u(3,i) - k(3)*u(2,i);
                    crosskv(2) = k(3)*u(1,i) - k(1)*u(3,i);
                    crosskv(3) = k(1)*u(2,i) - k(2)*u(1,i);
                    v_rotated(:,i) = cos(theta)*u(:,i) + (crosskv)*sin(theta)...
                        + k*(dot(k,u(:,i)))*(1 - cos(theta));
                end
            end
        end
        % GET THE UNIT VECTOR
        function [unitVector] = unit(vector)
            unitVector = vector/norm(vector);
        end
        % SKEW MATRIX
        function [V] = skew(v)
            assert(numel(v) == 3,'Input vector must be [3x1]');
            % SKEW THE VECTOR COMPONENTS
            V = [ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        end
    end
    %% EULER ROTATION MATHEMATICS
    methods (Static)
        % CONVERT BODY AXIS RATES TO EULER RATES
        function [eulerRates] = bodyAxisRatesToEularRates(eulerPosition,bodyAxisRates)
            % This function calculates the equivalent eular rotational
            % rates from defined body axis rates and a eular position.
            
            % INPUT HANDLING
            assert(numel(bodyAxisRates) == 3,'Incorrect number of body axis rates');
            assert(numel(eulerPosition) == 3,'Incorrect number of Euler rotations');
            
            % DEFINE THE CONVERSION MATRIX
            conversionMatrix = eye(3);
            conversionMatrix(1,2) = sin(eulerPosition(1))*tan(eulerPosition(2));
            conversionMatrix(1,3) = cos(eulerPosition(1))*tan(eulerPosition(2));
            conversionMatrix(2,2) = cos(eulerPosition(1));
            conversionMatrix(2,3) = -sin(eulerPosition(1));
            conversionMatrix(3,2) = sin(eulerPosition(1))*sec(eulerPosition(2));
            conversionMatrix(3,3) = cos(eulerPosition(1))*sec(eulerPosition(2)); % VALIDATED
            
            % DEFINE THE EULAR RATES FROM THE BODY AXIS RATES
            eulerRates = conversionMatrix*bodyAxisRates;
        end
        % CONVERT EULER RATES TO BODY AXIS RATES
        function [bodyAxisRates] = eularRatesToBodyAxisRates(eulerPosition,eulerRates)
            % This function calculates the conversion between the eular
            % rates and the body axis rates.
            
            % INPUT HANDLING
            assert(numel(eulerRates) == 3,'Incorrect number of Euler rates');
            assert(numel(eulerPosition) == 3,'Incorrect number of Euler rotations');
            
            % DEFINE THE CONVERSION MATRIX
            conversionMatrix = eye(3,3);
            conversionMatrix(1,3) = -sin(eulerPosition(2));
            conversionMatrix(2,2) =  cos(eulerPosition(1));
            conversionMatrix(2,3) =  cos(eulerPosition(2))*sin(eulerPosition(1));
            conversionMatrix(3,2) = -sin(eulerPosition(1));
            conversionMatrix(3,3) =  cos(eulerPosition(2))*cos(eulerPosition(1)); % VALIDATED
            % MAP THE EULAR ROTATION RATES TO THE BODY AXIS ROTATIONAL RATES
            bodyAxisRates = conversionMatrix*eulerRates;
        end
        % EULER ROTATION MATRIX TO EULER ANGLES
        function [eulers] = rotationMatrixToEulers(R_BG)
            % Rotation of the fixed to the global pose
            eulers = zeros(3,1);           
            eulers(1) =  atan2(R_BG(3,2),R_BG(3,3));  % The roll angle
            eulers(2) = -asin(R_BG(3,1));             % The pitch angle
            eulers(3) =  atan2(R_BG(2,1),R_BG(1,1));  % The yaw angle
        end
        % EULER ROTATION MATICES (STANDARD AEROSPACE 'ZYX' ROTATION ORDER)
        function [R] = eulersToRotationMatrix(eulers)
            % This rotation matrix maps the current rotated body axes to
            % the global axes.
            % http://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
            
            % This function converts the euler rotations (phi, theta, psi)
            % in the equivalent rotation matrix (R_AB) and inverse (R_BA).
            % "This matrix gives the conversion from body to navigation
            % frame"
            
            % ROLL AXIS ROTATIONS
            [R_phi] = OMAS_geometry.eulerToRotationMatrix_roll(eulers(1));            
            % PITCH AXIS ROTATIONS
            [R_theta] = OMAS_geometry.eulerToRotationMatrix_pitch(eulers(2));            
            % YAW AXIS ROTATIONS
            [R_psi] = OMAS_geometry.eulerToRotationMatrix_yaw(eulers(3));
            % BUILD THE 'ZYX' ROTATION MATRIX
            R = R_psi*R_theta*R_phi; % (Mapping local vectors to the global axes)
        end
        % EULER ROTATION MATRIX (Roll about the X-axis)
        function [R_phi] = eulerToRotationMatrix_roll(phi)
            % This rotation matrix maps a vector through the rotation angle
            % 'phi' about a reference 'x' axis.
            % ROLL AXIS ROTATIONS
            R_phi  = [ 1        0         0;
                       0 cos(phi) -sin(phi);
                       0 sin(phi)  cos(phi)];   
        end
        % EULER ROTATION MATRIX (Pitch about the Y-axis)
        function [R_theta] = eulerToRotationMatrix_pitch(theta)
            % This rotation matrix maps a vector through the rotation angle
            % 'theta' about a reference 'y' axis.
            % PITCH AXIS ROTATIONS
            R_theta = [cos(theta) 0 sin(theta);
                                0 1          0;
                      -sin(theta) 0 cos(theta)];
        end
        % EULER ROTATION MATRIX (Yaw about the Z-axis)
        function [R_psi] = eulerToRotationMatrix_yaw(psi)
            % This rotation matrix maps a vector through the rotation angle
            % 'psi' about a reference 'z' axis.
            % YAW AXIS ROTATIONS
            R_psi   = [cos(psi) -sin(psi) 0;
                       sin(psi)  cos(psi) 0;
                              0         0 1];
        end
        % GET THE ROTATION ANGLES BETWEEN TWO VECTORS
        function [angles,R]  = getVectorRotations(referenceVector,inputVector)
            % This function is designed to calculate the angles between a
            % given reference vector and an input vector; and a rotation
            % matrix between them.
            rotationAxis = cross(referenceVector,inputVector);
            rotationSkewMatrix = OMAS_geometry.skew(rotationAxis);
            s = abs(sqrt(sum(rotationAxis.^2)));         % Sin of angle
            c = dot(referenceVector,inputVector);        % Cos of angle
            R = eye(3) + rotationSkewMatrix + rotationSkewMatrix^2*((1-c)/s^2);
            % Get angles
            angles = abs([0;asin(s);acos(c)]);
        end
        % DEFINE THE PITCH AND YAW TO MEET TARGET HEADING
        function [lambda,theta] = getVectorHeadingAngles(V,U)
            % This function calculates the Line Of Sight (LOS) [horezontal]
            % angle which aids in defining the bank angle as a control input to the
            % dynamics.
            % INPUTS:
            % V - The current unity velocity vector
            % U - The unit correction vector
            % OUTPUTS:
            % lambda - The azimuth angle (2D)
            % theta  - The pitch angle   (3D)
                        
            % GET THE LINE OF SIGHT ANGLE             
            % NORMALISE THE ELEMENTS
            V = V/norm(V);
            U = U/norm(U);
            % GET THE HORIZONTAL ELEMENTS
            Vh = [V(1);V(2);0];
            Uh = [U(1);U(2);0];                                            % Reject the vertical elements
            % DEFINE LOS ANGLE            
            rotationAxis = cross(Vh,Uh);
            lambda = sign(rotationAxis(3))*acos(dot(Vh,Uh)/norm(Vh));      % Get the angle, signed by the direction of its cross product
            % GET THE ELEVATION ANGLE 
            theta = atan2(U(3),norm(Uh));
        end
    end
    %% QUATERNION ROTATION MATHEMATICS
    methods (Static)
        % ANALYTICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
        function [q] = getAnalyticalTriadRotation(referenceTriad,targetTriad)
            % This function defines the quaternion describing the rotations
            % between a reference triad and a second triad.
            
            % NORMALISE THE TRIAD VECTORS
            for dim = 1:size(referenceTriad,2)
                [referenceTriad(:,dim)] = OMAS_geometry.unit(referenceTriad(:,dim));
                [targetTriad(:,dim)] = OMAS_geometry.unit(targetTriad(:,dim));
            end
            
            % EXTRACT REFERENCE AXES
            xAxis = targetTriad(:,1);
            zAxis = targetTriad(:,3);   % Get the rotated body axes
            xAxis_ref = referenceTriad(:,1);
            zAxis_ref = referenceTriad(:,3); % Get the reference unit ENU triad (trying to get to)

            % FIRST ALIGN THE Z-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GET THE QUATERNION ROTATION TO ALLIGN THE Z-AXES
            [q_zAlign] = OMAS_geometry.vectorsToQuaternion(zAxis,zAxis_ref);    % Quaternion aligning global and body z vectors
            [R_zAlign] = OMAS_geometry.quaternionToRotationMatrix(q_zAlign);    % Equivalent rotation matrix
            % ALIGN THE X AXIS IN THE Z-PLANE
            xAxis_intermediate = R_zAlign*xAxis;
            % TAKE ITS PROJECTIONS IN THE XY PLANE & RENORMALISE
            xAxis_intermediate(3) = 0;
            xAxis_intermediate = OMAS_geometry.unit(xAxis_intermediate);
            % OTHERWISE JUST ALIGN THE X-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GET THE QUATERNION ROTATION TO ALLIGN THE X-AXES
            [q_xAlign] = OMAS_geometry.vectorsToQuaternion(xAxis_intermediate,xAxis_ref);
            [R_xAlign] = OMAS_geometry.quaternionToRotationMatrix(q_xAlign);
            
            % COMPUTE THE COMPOSITE ROTATION MATRIX
            comp_rotation = R_xAlign * R_zAlign;
            % COVERT THE ROTATION MATRIX TO A QUATERNION DESCRIBING THE
            % ROTATION FROM TRIAD_REF TO TRIAD_FINAL
            q = rotm2quat(comp_rotation)';
%             q = OMAS_geometry.rotationMatrixToQuaternion(comp_rotation);
            % OUTPUT CATCHA
            assert(any(isnan(q)) == 0,'No valid quaternion found');
        end
        % NUMERICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
        function [q] = getNumericalTriadRotation(referenceTriad,targetTriad)
            % This function defines the quaternion nescessary to rotate one
            % axis triad into alignment with a second triad. This is done by
            % solving the rotation matrix-defined symbolic expressions for
            % the equivalent quaternion rotation.
            
            syms q1 q2 q3 q4 real
            
            % Define the symbolic quaternion rotation matrix
            R_q = zeros(3);
            R_q(1,1) = q1^2 + q2^2 - q3^2 - q4^2;
            R_q(1,2) = 2*(q1*q4 + q2*q3);
            R_q(1,3) = 2*(q2*q4 - q1*q3);
            R_q(2,1) = 2*(q3*q2 - q1*q4);
            R_q(2,2) = q1^2 - q2^2 + q3^2 - q4^2;
            R_q(2,3) = 2*(q1*q2 + q3*q4);
            R_q(3,1) = 2*(q1*q3 + q2*q4);
            R_q(3,2) = 2*(q3*q4 - q1*q2);
            R_q(3,3) = q1^2 - q2^2 - q3^2 + q4^2;
            % DEFINE THE SYMBOLIC EXPRESSIONS
            expressions = [R_q*referenceTriad(:,1) - targetTriad(:,1);
                R_q*referenceTriad(:,2) - targetTriad(:,2);
                R_q*referenceTriad(:,3) - targetTriad(:,3);
                q1^2 + q2^2 + q3^2 + q4^2 - 1];
            expressions = symfun(expressions,[q1 q2 q3 q4]);                % Define symbolic function in q_{1-4}
            modeqn = @(x)double(expressions(x(1),x(2),x(3),x(4)));          % Rearrange to be in terms of q1,q2 etc
            % DEFINE SOLVER OPTIONS
            OPTIONS = optimoptions('fsolve','Algorithm','levenberg-marquardt');
            % PASS SYMBOLIC EQUATIONS TO F-SOLVE
            x0 = [ 1 0 0 0];                 % INITIAL QUATERNION GUESS
            q = fsolve(modeqn,x0,OPTIONS);   % SOLVE FOR THE ROTATION QUATERNION
            % RE-NORMALISE THE QUATERNION
            [q] = OMAS_geometry.qUnit(q);
        end
        % UPDATE A GIVEN QUATERNION
        function [q_new] = integrateQuaternion(q0,omega,dt)
            % This function is designed to compute the quaternion update routine from
            % one quaternion angle to all the appropriate representations.
            % INPUTS
            % q         - The initial quaternion axis reference
            % omega - The current axis rates [p;q;r]
            % dt        - The timestep between instances
            % OUTPUTS
            % q         - The updated quaternion
            % q_dot     - The quaternion difference
            % R_q       - The rotation matrix R_G to R_B
            
            % GET THE QUATERNION DIFFERENTIAL
            [q_dot] = OMAS_geometry.quaternionDifferential(q0,omega);    
                      
            % Integrate the quaternion difference
            q_new = q_dot*dt + q0;
            % Re-normalise
            q_new = OMAS_geometry.qUnit(q_new);
        end             % NOTATION B
        % GET THE QUATERNION DIFFERNTIAL (as a function of omega)
        function [q_dot] = quaternionDifferential(q0,omega)
            % This function was created to allow the quaternion
            % differential to be called seperately to in the integration
            % function. This allows the quaternion 'q0' to be integrated
            % externally.
            
            % Calculate the integration drift normalising elements
            K = 1; % k*dt <= 1                                         
            lambda = 1 - (q0(1)^2 + q0(2)^2 + q0(3)^2 + q0(4)^2);
         
            % Calculate the quaternion difference                          % NOTATION A
%             q_dot = 0.5*[K*lambda,  omega(3), -omega(2), omega(1);
%                     -omega(3),      K*lambda,  omega(1), omega(2);
%                      omega(2), -omega(1),      K*lambda, omega(3);
%                     -omega(1), -omega(2), -omega(3),     K*lambda]*q0;

            % Calculate the quaternion difference                          % NOTATION B
            q_dot = 0.5*[ K*lambda, -omega(1), -omega(2), -omega(3);
                          omega(1),  K*lambda,  omega(3), -omega(2);
                          omega(2), -omega(3),  K*lambda,  omega(1);
                          omega(3),  omega(2), -omega(1),  K*lambda]*q0;                
        end
        % ROTATE A VECTOR THROUGH A QUATERNION
        function [newVector] = quaternionVectorRotation(q,oldVector)
            % This function rotates a cartesian vector through a quaternion
            % rotation. Associated block:
            % "Quaternion Rotation"
            if numel(q) ~= 4
                %                 error('Input quaternion is of length %s',num2str(length(q)));
                error(char(['Input quaternion is of length ',num2str(length(q))']));
            end
            % Normalise the quaternion rotation
            q = OMAS_geometry.qUnit(q);
            % Rotate the vector through the quaternion elements
            newVector = 2*[(0.5 - q(3)^2 - q(4)^2), (q(1)*q(4) + q(2)*q(3)), (q(2)*q(4) - q(1)*q(3));...
                (q(2)*q(3) - q(1)*q(4)), (0.5 - q(2)^2 - q(4)^2), (q(1)*q(2) + q(3)*q(4));...
                (q(1)*q(3) + q(2)*q(4)), (q(3)*q(4) - q(1)*q(2)), (0.5 - q(2)^2 - q(3)^2)]*oldVector;
        end
        % GET THE QUATERNION BETWEEN TWO VECTORS
        function [q] = vectorsToQuaternion(u,v)
            q = zeros(4,1);
            % Normalise the quaternion
            u = u/norm(u);
            v = v/norm(v);
            % Get the axis vector
            q(2:4) = cross(u,v);
            % Define the rotation about that vector
            q(1) = sqrt((norm(u)^2)*(norm(v)^2)) + dot(u,v);
            % Normalise the quaternion
            [q] = OMAS_geometry.qUnit(q);
        end
        % CONVERT ROTATION ANGLES INTO A QUATERNION
        function [q] = eulersToQuaternion(eulerAngles)
            % THis function converts euler angle rotations into a
            % quaternion via its components. Associated block:
            % "Rotation Angles to Quaternions".
            eulerAngles = 0.5*eulerAngles;
            % Build vector of trigonometric arguements
            trigArgs = [sin(eulerAngles);cos(eulerAngles)];
            % Assemble Quaternion components
            q = zeros(4,1);
            q(1) = trigArgs(4)*trigArgs(5)*trigArgs(6) + trigArgs(1)*trigArgs(2)*trigArgs(3);
            q(2) = trigArgs(4)*trigArgs(5)*trigArgs(3) - trigArgs(1)*trigArgs(2)*trigArgs(6);
            q(3) = trigArgs(4)*trigArgs(2)*trigArgs(6) + trigArgs(1)*trigArgs(5)*trigArgs(3);
            q(4) = trigArgs(1)*trigArgs(5)*trigArgs(6) - trigArgs(4)*trigArgs(2)*trigArgs(3);
        end
        % CONVERT QUATERNION INTO ROTATION ANGLES
        function [eulerAngles] = quaternionToEulers(q)                     % NOTATION B
            % Gets the rotation angles from an equivalent quaternion
            % attitude. Associated block:
            % "Quaternions to Rotation Angles"
            if length(q) ~= 4
                error('Input quaternion is of length %s',num2str(length(q)));
            end
            % Normalise the quaternion
            q = OMAS_geometry.qUnit(q);
            eulerAngles = zeros(3,1);
            % QUATERNION NOTATION B [roll pitch yaw]
            eulerAngles(1) = atan2(2*(q(1)*q(2) + q(3)*q(4)),(1 - 2*(q(2)^2 + q(3)^2)));              % Phi   (roll)
            eulerAngles(2) = asin( 2*(q(1)*q(3) - q(4)*q(2)));             % Theta (pitch)              
            eulerAngles(3) = atan2(2*(q(1)*q(4) + q(2)*q(3)),(1 - 2*(q(3)^2 + q(4)^2)));              % Psi   (yaw)
        end
        % CONVERT ROTATION MATRIX TO QUATERNION
        function [q] = rotationMatrixToQuaternion(R)
            % This function is designed to convert from a rotation matrix
            % to an equivalent quaternion. This function is also parallel
            % to "rotm2quat.m".
            
            q = zeros(4,1);
            % Assemble the quaternion elements
            q(1) = 0.5*sqrt(1 + R(1,1) + R(2,2) + R(3,3));
            q(2) = (R(3,2) - R(2,3))/(4*q(1));
            q(3) = (R(1,3) - R(3,1))/(4*q(1));
            q(4) = (R(2,1) - R(1,2))/(4*q(1));
            
            assert(any(isnan(q)) == 0,'quaternion calculation failed');

            % Normalise the quaternion
            q = OMAS_geometry.qUnit(q);
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R] = quaternionToRotationMatrix(q)
            % This function defines the equivalent rotation matrix from the
            % provided quaternion. (Associated block "Quaternion to DCM")
            % INPUT:
            % q    - The quaternion rotation
            % OUTPUT:
            % R_AB - The rotation matrix through A to B
            % R_BA - The reverse rotation matrix from B to A
            
            assert(numel(q) == 4,'Input quaternion has wrong dimensions');
            
            R = zeros(3,3);
            % Normalise the quaternion
            q = OMAS_geometry.qUnit(q);                                   % NOTATION A
            % Assemble the quaternion elements
%             R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
%             R(1,2) = 2*(q(1)*q(4) + q(2)*q(3));
%             R(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
%             R(2,1) = 2*(q(3)*q(2) - q(1)*q(4));
%             R(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
%             R(2,3) = 2*(q(1)*q(2) + q(3)*q(4));
%             R(3,1) = 2*(q(1)*q(3) + q(2)*q(4));
%             R(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
%             R(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
            
            % SAME AS THE ROBOTICS TOOL BOX....
            R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;                 % NOTATION B
            R(1,2) = 2*(q(2)*q(3) - q(1)*q(4));
            R(1,3) = 2*(q(1)*q(3) + q(2)*q(4));
            R(2,1) = 2*(q(2)*q(3) + q(1)*q(4));
            R(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
            R(2,3) = 2*(q(3)*q(4) - q(1)*q(2));
            R(3,1) = 2*(q(2)*q(4) - q(1)*q(3));
            R(3,2) = 2*(q(1)*q(2) + q(3)*q(4));
            R(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
        end
        
        % THE QUATERNION DIFFERENCE
        function [dq] = qDifference(q1,q2)
            % This function gets the quaternion that will rotate from the
            % q1 orientation to the q2 orientation [q1->q2].
            
            assert(size(q1,1) == 4,'The initial quaternion (q1) must be a 4x1 column vector')
            assert(size(q1,1) == 4,'The terminal quaternion (q2) must be a 4x1 column vector')
            
            % Get the quaternion inverse
            [q1_inv] = OMAS_geometry.qInverse(q1);
            % Calculate difference
            dq = OMAS_geometry.qMultiply(q1_inv,q2);
            % Re-normalise
            dq = OMAS_geometry.qUnit(dq);
        end       
        % QUATERNION DIVISION
        function [qDiv] = qDivide(q,r)
            % The function compute the division of quaternion q, by
            % quaterion v.
            % Associated block:
            % "Quaternion Division"
            
            qDiv = zeros(4,1);
            % Normalise the second vector
            qDiv(1) =  q(1)*r(1) + q(2)*r(2) + q(3)*r(3) + q(4)*r(4);
            qDiv(2) = -q(1)*r(2) + q(2)*r(1) + q(3)*r(4) - q(4)*r(3);
            qDiv(3) =  q(1)*r(3) + q(2)*r(4) + q(3)*r(1) + q(4)*r(2);
            qDiv(4) = -q(1)*r(4) + q(2)*r(3) - q(3)*r(2) + q(4)*r(1);
            qDiv = qDiv./(OMAS_geometry.qNorm(r));
        end
        % QUATERNION MULTIPLICATION
        function [qv] = qMultiply(q,v)
            % Calculate the product of two quaternions
            % Associated block:
            % "Quaternion Multiplication"
            % Multiply the quaternion elements
            
            assert(size(q,1) == 4 && size(v,1) == 4,...
                'Both quaternion must be provided as 4x1 column vectors')
            % Quaternion projection matrix
            qv = [v(1), -v(2), -v(3), -v(4);
                  v(2),  v(1), -v(4),  v(3);
                  v(3),  v(4),  v(1), -v(2);
                  v(4), -v(3),  v(2),  v(1)]*q; % Confirmed with matlab website
            % Re-normalise the output
            qv = OMAS_geometry.qUnit(qv);
        end
        % QUATERNION INVERSE
        function [qInv] = qInverse(q)
            % The quaternion norm
            qSqr = q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2;
            % CALCULATE THE INVERSE OF A GIVEN QUATERNION
            qInv = zeros(4,1);
            qInv(1) =  q(1)/qSqr;
            qInv(2) = -q(2)/qSqr;
            qInv(3) = -q(3)/qSqr;
            qInv(4) = -q(4)/qSqr;
        end
        % QUATERNION CONJUGATE
        function [qCon] = qConjugate(q)
            % This function defines the conjugate of the given quaternion.
            % Associated block:
            % "Quaternion Conjugate"
            qCon = zeros(4,1);
            qCon(2) = -q(2);
            qCon(3) = -q(3);
            qCon(4) = -q(4);
        end
        % THE UNIT QUATERNION
        function [qUnit] = qUnit(q)
            qUnit = q./OMAS_geometry.qNorm(q);
        end
        % THE QUATERNION NORM
        function [qNorm] = qNorm(q)
            % GETS THE NORM OF THE QUATERNION VECTOR
            qNorm = sqrt(sum(q.^2));    % VALIDATED
        end
    end
    
    %% 3D GEOMETRIC TOOLS
    methods (Static)
        % INTERSECT TWO 3D CUBOIDS
        function [intersectFlag] = cuboidIntersect(minA,maxA,minB,maxB)
            % This function determines if two 3D cuboids currently intersect
            for dim = 1:numel(minA)
                % FOR A COLLISION, THE CUBOIDS MUST INTERSECT IN ALL DIMENSIONS
                if ~OMAS_geometry.rangeIntersect(minA(dim),maxA(dim),minB(dim),maxB(dim))
                    intersectFlag = 0;
                    return
                else
                    intersectFlag = 1;
                end
            end
        end
        % CHECK IF A POINT LIES WITHIN A RANGE
        function [flag] = rangeIntersect(minA,maxA,minB,maxB)
            assert(minA < maxA,'Minimum value must be less than the maximum');
            assert(minB < maxB,'Minimum value must be less than the maximum');

            flag = 0;
            if maxA >= minB && minA <= maxB
                flag = 1;
            end        
        end
        
        % THE EXPRESSOIN FOR A PLANE DEFINED BY POINTS
        function [a,b,c,d] = planarExpression(p1,p2,p3)
            % This function defines the parameteric equation of a plane
            % from three defining vertex points. The result is of the form:
            % ax + by + cz + d = 0;

            a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) ...
                - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) );
            b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) ...
                - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) );
            c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) ...
                - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) );
            d = - p2(1) * a - p2(2) * b - p2(3) * c;
        end
    end
    
    %% UNIVERSAL DRAWING MECHANISM
    methods (Static)
        % DRAW SPHERE
        function [sphereHandle] = drawSphere(figureHandle,position,radius)
            % This function returns a sphere coordinate cloud at a given
            % position in space, of a given radius. 'figureHandle' is used
            % to provide the function context.
            [X,Y,Z] = sphere(40);
            X = X.*radius + position(1);
            Y = Y.*radius + position(2);
            Z = Z.*radius + position(3);
            % DEFINE THE SPHERE MESH
            sphereHandle = mesh(gca,X,Y,Z);
            set(sphereHandle,...
                'FaceAlpha',0.25,...
                'LineWidth',0.1,...
                'EdgeAlpha',0.1);
        end
        % DRAW UNIT TRIAD
        function [triadVectors] = drawTriad(figureHandle,position,R)
            % Draw a unit triad at a cartesian position, rotated by R.
            % FigureHandl Here is used to bring the current figure to the
            % function context.
            
            colourVector = {'r','g','b'};
            triadVectors = eye(3);
            for axis = 1:size(triadVectors,2)
                triadVectors(:,axis) = R*triadVectors(:,axis);
                % Draw vectors
                quiver3(gca,position(1),position(2),position(3),triadVectors(1,axis),triadVectors(2,axis),triadVectors(3,axis),colourVector{axis});
            end
        end 
    end
end

