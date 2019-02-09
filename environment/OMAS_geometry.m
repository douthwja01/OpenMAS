%% OPENMAS 3D-AXIS MANIPULATION TOOLS (OMAS_geometry.m) %%%%%%%%%%%%%%%%%%%
% This class provides the axis manipulation utilities used by OpenMAS to
% evaluate object progression over time.

% Author: James A. Douthwaite 06/09/2017

classdef OMAS_geometry
    % This class is designed to handle all the axis transformations and
    % conversions as part of the OpenMAS simulator.
    
    %% GENERIC TOOLS
    methods (Static) 
        % PLANAR PROJECTION OF A VECTOR ON A PLANE
        function [v_proj] = vectorPlanarProjection(n,v)
            % Normalise the normal vector
            n = OMAS_geometry.unit(n);
            % Get the projection on the plane
            v_proj = v - (dot(n,v)/norm(n))*n; 
        end
        % SIGNED ANGLE BETWEEN TO VECTORS
        function [planarAngle,isRightLogical] = signedPlanarAngle(n,u,v)
            % This function computes the signed angle between two vectors
            % with respect to a provided normal.
            
            % RENORMALISE THE INPUTS
            n = OMAS_geometry.unit(n);
            u = OMAS_geometry.unit(u);
            v = OMAS_geometry.unit(v);
            % DETERMINE THE DIRECTION (CW+ve with respect to the normal)
            isRightLogical = -sign(dot(cross(u,v),n));
            % DETERMINE THE MAGNITUDE
            [planarAngle] = OMAS_geometry.planarAngle(u,v);
            % COMBINE
            planarAngle = isRightLogical*planarAngle;
        end
        % PLANAR ANGLE BETWEEN TWO VECTORS
        function [planarAngle] = planarAngle(u,v)
            planarAngle = acos(dot(u,v)/(norm(u)*norm(v)));
        end
        % ROTATE A VECTOR ABOUT AN AXIS VECTOR (RODRIGUES)
        function [v_rotated] = rodriguesRotation(u,k,theta)
            % v - Vector to be rotated
            % k - Is the rotation axis
            % Theta - The angle the vector is to be rotated through
            
            assert(numel(u) == 3,'Rotation vector must be of size [3x1].')
            assert(numel(k) == 3,'The rotation axis must be of size [3x1]');
            assert(numel(theta) == 1,'The rotation angle %.0f must be a scalar',theta);
            
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
        % NORMALISE ANGLE BETWEEN PI & -PI
        function [angle] = normaliseAngle(angle)
           % This function normalise an angle between -pi and pi.
           
           piIntegers = angle/2*pi;
           
           if piIntegers > 1 
               angle = angle - piIntegers*(2*pi);
           elseif piIntegers < -1
               angle = angel + piIntegers*(2*pi);
           end
%             if angle > pi
%                 angle = angle - 2*pi;
%             end
%             if angle < -pi
%                 angle = angle + 2*pi;
%             end
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
            
            % RENORMALISE
            q = OMAS_geometry.qUnit(q);
            % OUTPUT CATCHA
            assert(any(isnan(q)) == 0,'No valid quaternion found');
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
                error(char(['Input quaternion is of length %d ',length(q)']));
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
            assert(size(q,1) == 4,'A quaternion column vector is expected 4x1');

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
        % INTERSECT - SPHERE AND ROTATED CUBOID (OBB)
        function [intersectFlag,separation_OBB]  = intersect_OBB_sphereCuboid(centerA,radiusA,centerB,cuboidVerticesB)
            % This function computes the intersection between a sphere and a 
            % rotated bounding box (OBB). This function assumes that the
            % problem is being resolved in the axes of A (i.e) 'R' is the
            % rotation matrix taking the cube from its own frame to the
            % axes of A.
            
            % RELATIVE POSITION RAY
            [BAray] = OMAS_geometry.defineRay(centerB,(centerA-centerB));       % Vector between             
            % If we define A to be sphere, its projection on the separation
            % vector will always be its radius
            projection_A = radiusA;            
            % GET THE PROJECTION OF ALL THE VERTICES ON THE SEPARATING VECTOR
            projection_B = zeros(size(cuboidVerticesB,1),1);
            for i = 1:size(cuboidVerticesB,1)
                % Its projection on the separation vector
                [projection_B(i)] = OMAS_geometry.pointProjectionOnRay(BAray,cuboidVerticesB(i,:)');
            end
            % THE MAXIMAL PROJECTION TOWARDS THE SPHERE
            projection_B = max(projection_B);
            % CHECK IF THIS EXCEEDS THE SEPARATION OF THE TWO OBJECTS
            separation_OBB = norm(centerB-centerA) - (projection_A + projection_B);      
            % EVALUATE INTERSECTION
            intersectFlag = 0;
            if separation_OBB < 0
                separation_OBB = 0;
                intersectFlag = 1;
            end
        end
        % INTERSECT - TWO ROTATED CUBOIDS (OBB)
        function [intersectFlag,separation_OBB]  = intersect_OBB_cuboids(centerA,cuboidVerticesA,centerB,cuboidVerticesB)
            
            % RELATIVE POSITION RAYS
            ABdirection = (centerB-centerA)/norm(centerB-centerA);
            [ABray] = OMAS_geometry.defineRay(centerA, ABdirection);
            [BAray] = OMAS_geometry.defineRay(centerB,-ABdirection);
            
            % PROJECTIONS OF EACH VERTEX ONTO THE SEPARATION VECTOR
            projections_A = zeros(size(cuboidVerticesA,1),1);
            for i = 1:size(cuboidVerticesA,1)
                % Its projection on the separation vector
                [projections_A(i)] = OMAS_geometry.pointProjectionOnRay(ABray,cuboidVerticesA(i,:)');
            end
            projections_B = zeros(size(cuboidVerticesB,1),1);
            for i = 1:size(cuboidVerticesB,1)
                % Its projection on the separation vector
                [projections_B(i)] = OMAS_geometry.pointProjectionOnRay(BAray,cuboidVerticesB(i,:)');
            end
            % MAXIMAL PROJECTIONS ALONG THE SEPARATION VECTOR
            tmax_A = max(projections_A);
            tmax_B = max(projections_B);
            % CHECK IF THIS EXCEEDS THE SEPARATION OF THE TWO OBJECTS
            separation_OBB = norm(centerB-centerA) - (tmax_A + tmax_B);
            % EVALUATE INTERSECTION
            intersectFlag = 0;
            if separation_OBB < 0
                separation_OBB = 0;                                            % Catch negative separations
                intersectFlag = 1;
            end
        end
        % INTERSECT - SPHERE AND AXIS ALIGNED CUBOID (AABB)
        function [intersectFlag,separation_AABB] = intersect_AABB_sphereCuboid(center,radius,cubeMins,cubeMaxs)
            % get cuboid closest point to sphere center by clamping
            %             x = max(cubeMins(1), min(sphere.x, cubeMaxs(1)));
            %             y = max(cubeMins(2), min(sphere.y, cubeMaxs(2)));
            %             z = max(cubeMins(3), min(sphere.z, cubeMaxs(3)));
            
            x = max([cubeMins(1), min([center(1), cubeMaxs(1)])]);
            y = max([cubeMins(2), min([center(2), cubeMaxs(2)])]);
            z = max([cubeMins(3), min([center(3), cubeMaxs(3)])]);
            
            % this is the same as isPointInsideSphere
            separation_AABB = sqrt((x - center(1)) * (x - center(1)) + ...
                                   (y - center(2)) * (y - center(2)) + ...
                                   (z - center(3)) * (z - center(3)));
            % IS CUBOID VERTEX (CLOSEST TO SPHERE CENTER) INSIDE SPHERE
            intersectFlag = separation_AABB < radius;
        end
        % INTERSECT - AXIS ALIGNED CUBOIDS (AABB)
        function [intersectFlag] = intersect_AABB_cuboids(minA,maxA,minB,maxB)
            % This function determines if two 3D cuboids currently intersect
            intersectFlag = 1;
            iMatrix = zeros(numel(minA),2);
            for dim = 1:numel(minA)
                % INTERSECTING TWO RANGES, EACH DIMENSION AT A TIME
                [iMatrix(dim,1),iMatrix(dim,2)] = OMAS_geometry.rangeIntersect(minA(dim),maxA(dim),minB(dim),maxB(dim));
                if ~iMatrix(dim,1)
                    intersectFlag = 0;
                end
            end
            % THE MEAN SEPARATION
            %separation_AABB = sqrt(sum(iMatrix(:,2).^2));
        end
        % INTERSECT - SPHERES
        function [intersectFlag] = intersect_spheres(positionA,radiusA,positionB,radiusB)
            % This uses a simple one dimensional comparison to determine
            % the overlap of two spherical volumes.
            
            % COLLISION LOGIC
            intersectFlag = 0;
            if norm(positionB - positionA) - (radiusA + radiusB) < 0
                intersectFlag = 1;
            end
        end
        % INTERSECT - RANGE
        function [intersectFlag,separation_range]  = rangeIntersect(minA,maxA,minB,maxB)
            % This function determines the intersection between two 1D
            % ranges, and the distance between them.
            assert(minA < maxA,'Minimum value must be less than the maximum');
            assert(minB < maxB,'Minimum value must be less than the maximum');
            % Do either of the ranges on an assumed axis overlap
            if maxA <= minB 
                intersectFlag = 0;              % They intersect
                separation_range = minB - maxA;
            elseif minA >= maxB
                intersectFlag = 0;              % They intersect
                separation_range = minA - maxB;
            else
                intersectFlag = 1;              % They don't
                separation_range = 0;
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
        % THE SQUARED DISTNACE TO A LINE SEQMENT
        function [dSq] = distSqPointLineSegment(p1,p2,q)
            % This function computes the squared distnace from a line
            % segment with the specified endpoints to a specified point.
            % INPUTS:
            % p1 - The first endpoint of a line segment
            % p2 - The second endpoint of a line segment
            % q  - The point from which the squared distance is measured
            
            % Compute 
            r = ((q - p1) * (p2 - p1)) / absSq(p2 - p1);
            
            if (r < 0) 
                dSq = (q - p1)^2;
            elseif (r > 1) 
                dSq = (q - p2)^2;
            else
                dSq = (q - (p1 + r * (p2 - p1)))^2;
            end
            
        end
        % COMPUTE THE SIGNED DISTANCE FROM A LINE CONNECTING POINTS TO POINT
        function [d_signed] = leftOf(a,b,c)
           % This function computes the signed distance from a line 
           % connecting the specified points to a specified point.
           % INPUTS:
           % p1 - The first point on the line
           % p2 - The second point on the line
           % q  - The point from which the signed distance is measured
           % OUTPUT:
           % d_signed - Positive when the point p3 lies to the left of the
           %            line p1-p2
           
           % The signed distance betweeen the points
           d_signed =  OMAS_geometry.det((a - c),(b - a));
        end
        % CALCULATE THE DETERMINANT 
        function [d] = det(p1,p2)
            % Calculate the determinant of the dimensional square matrix.
            d = p1(1)*p2(2) - p1(2)*p2(1);
        end
        
        % /////////////////////// RAY FUNCTIONS ///////////////////////////
        % RAY-SPHERE INTERSECTION
        function [p_int] = raySphereIntersection(ray,centroid,radius)
           % This function calculates the point(s) of intersection between 
           % a ray and a defined sphere. 
           
           % INPUTS:
           % ray      - The 
           % centroid - The center of the defined sphere.
           % radius   - The radius of the defined sphere.
           
           % OUTPUT CONTAINER
           p_int = [];
           tolerance = (1E-3);
           
           % The projection of the separation on the ray direction
           t = dot((centroid-ray.origin),ray.direction);
           
           if t < 0 
               % Ray is directed away from the sphere
               return
           end
           
%            if t > ray.magnitude
%                % The ray does not extend that far
%                return
%            end
           
           % Rays minimal radius
           ySq = t^2 - (norm(centroid-ray.origin))^2;
           
           % Get distance between the centerline and the point on the
           % circumference
           rSq = radius^2;
           
           % If the radial distance greater than the projection on the radius
           if ySq < rSq
               % The magnitude difference
               dt = sqrt(rSq + ySq);
               % First intersection distance
               t1 = t - dt;
               % Require a minimum segment length
               if dt < tolerance 
                  p_int = (ray.origin + ray.direction*t1)';
                  return
               end
               % Second intersection distance
               t2 = t + dt;
               % The second
               p_int = vertcat(p_int,(ray.origin + ray.direction*t2)');
           end
        end
        % PROJECTION OF A RAY ON A PLANE DEFINED BY A NORMAL
        function [v_prj] = rayProjectionOnPlane(ray,planeNormal)
            % This function computes the projection of a vector on a plane
            % defined by its normal vector.
            
            % Ensure normal vector is the normal
            [planeNormal] = OMAS_geometry.unit(planeNormal);
            
            if isfield(ray,'magnitude')
                % Scale the ray to its maximum length
                v_prj = ray.magnitude*ray.direction - dot(ray.direction,planeNormal)*planeNormal; 
            else
                % The projection of the ray on the plane
                v_prj = ray.direction - dot(ray.direction,planeNormal)*planeNormal;
            end
            
        end
        % PROJECTION OF A POINT ON A VECTOR
        function [v_mag] = pointProjectionOnRay(ray,point)
            % Calcuate the projection magnitude(distance)
            v_mag = dot(ray.direction,(point - ray.origin));
        end
        % DEFINE RAY
        function [ray] = defineRay(origin,direction,magnitude)
           % Define a ray with a given a origin and direction.
           assert(numel(origin) == 3 && numel(direction) == 3,'Parameters must be 3D');
           if nargin < 3
               magnitude = 1;
           end
           % NORMALISE THE DIRECTION VECTOR
           direction = direction/norm(direction);
           % DEFINE THE RAY STRUCTURE
           ray = struct('origin',origin,...
                        'direction',direction,...
                        'magnitude',magnitude);
        end
        % PLOT RAY
        function plotRay(ray)
            % This function simply plots a representation of the ray to the
            % current axes for visualisation
            
            % INPUT CHECKING
            assert(numel(ray.direction) == numel(ray.origin),'The dimensions of the ray are inconsitance');
            
            if isfield(ray,'magnitude')
                direction = ray.direction*ray.magnitude;
            else
                direction = ray.direction;
            end
            
            colour = 'k';
            
            % PLOT BASED ON THE DIMENSIONALITY
            if numel(ray.direction) == 3
                % THE RAY ORIGIN
                scatter3(gca,ray.origin(1),ray.origin(2),ray.origin(3),[colour 'o']);
                % THE RAY DIRECTION
                q = quiver3(gca,ray.origin(1),ray.origin(2),ray.origin(3),...
                            direction(1),direction(2),direction(3),colour);
            else
                % THE RAY ORIGIN
                scatter(gca,ray.origin(1),ray.origin(2),[colour 'o']);
                % THE RAY DIRECTION
                q = quiver(gca,ray.origin(1),ray.origin(2),...
                           direction(1),direction(2),colour);
            end
            q.AutoScaleFactor = 1;
        end
    end
end

