%% OPENMAS AXIS MANIPULATION SUITE (OMAS_axisTools.m) %%%%%%%%%%%%%%%%%%%%%
% This class provides the axis manipulation utilities used by OpenMAS to
% evaluate object progression over time.

% Author: James A. Douthwaite 06/09/2017

classdef OMAS_axisTools
    % This class is designed to handle all the axis transformations and
    % conversions as part of the OpenMAS simulator.
    
    methods (Static) % METHODS WITHOUT SELF REFERENCE   
        %% GENERIC TOOLS
        % CONVERT NED 12DOF STATE TO ENU
        function [FLUstate] = convertNEDstateToENU(FRDstate)
            % This function simply converts an NED state vector into its
            % equivalent (aligned) ENU state vector.
            
            FLUstate = zeros(12,1);
            % POSITIONS
            FLUstate(1) =  FRDstate(1);   % Forward motion:   + East motion
            FLUstate(2) = -FRDstate(2);   % Rightward motion: - Northward motion
            FLUstate(3) = -FRDstate(3);   % Downward motion:  - Upward motion
            % ROTATIONS
            FLUstate(4) =  FRDstate(4);    
            FLUstate(5) =  FRDstate(5);   
            FLUstate(6) = -FRDstate(6);  % Axis rotation is opposite    (looking along the axis)
            % LINEAR VELOCITIES
            FLUstate(7) =  FRDstate(7);
            FLUstate(8) = -FRDstate(8);   
            FLUstate(9) = -FRDstate(9);   % Same as positions
            % ANGULAR VELOCITIES
            FLUstate(10) =  FRDstate(10);
            FLUstate(11) =  FRDstate(11);  
            FLUstate(12) = -FRDstate(12); % Same as rotations
        end
        % CONVERT NED 12DOF STATE TO MATLAB (XYZ)
        function [XYZstate] = convertNEDstateToMatlab(NEDstate)
            XYZstate = zeros(12,1);
            % POSITIONS
            XYZstate(1) = NEDstate(1); % +x = +East
            XYZstate(2) = NEDstate(2); % +y = +North 
            XYZstate(3) = NEDstate(3); % +z = +Up
            % ROTATIONS
            XYZstate(4) = NEDstate(4);   % Rotations in matlab -ve  
            XYZstate(5) = NEDstate(5);    % X/Y axis directions inverted
            XYZstate(6) = NEDstate(6);   % Rotations in matlab -ve 
            % LINEAR VELOCITIES
            XYZstate(7) = NEDstate(7);
            XYZstate(8) = NEDstate(8);   
            XYZstate(9) = NEDstate(9); % As positions
            % ANGULAR VELOCITIES
            XYZstate(10) = NEDstate(10);
            XYZstate(11) = NEDstate(11);
            XYZstate(12) = NEDstate(12);
        end
        % CONVERT MATLAB XYZ 12DOF STATE TO ENU
        function [XYZstate] = convertENUStateToMatlab(ENUstate)
            % This function simply converts an NED state vector into its
            % equivalent (aligned) ENU state vector.
            
            XYZstate = zeros(12,1);
            % POSITIONS
            XYZstate(1) = ENUstate(1); % +x = +East
            XYZstate(2) = ENUstate(2); % +y = +North 
            XYZstate(3) = ENUstate(3); % +z = +Up
            % ROTATIONS
            XYZstate(4) = -ENUstate(5);   % Rotations in matlab -ve  
            XYZstate(5) = ENUstate(4);    % X/Y axis directions inverted
            XYZstate(6) = -ENUstate(6);   % Rotations in matlab -ve 
            % LINEAR VELOCITIES
            XYZstate(7) = ENUstate(7);
            XYZstate(8) = ENUstate(8);   
            XYZstate(9) = ENUstate(9); % As positions
            % ANGULAR VELOCITIES
            XYZstate(10) = -ENUstate(11);
            XYZstate(11) = ENUstate(10);
            XYZstate(12) = -ENUstate(12);
        end
        
        %% GENERIC AXIS TOOLS
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
        % CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
        function [outputMatrix] = skew(inputVector)
            % This function generates a skew-symmetric for the computation
            % of the vector cross-product.
            % INPUTS:
            % inputVector - The original 3D vector
            % OUTPUT:
            % outputMatrix - The equivalent skew-symmetric matrix
            
            if length(inputVector) ~= 3
%                 warning('The input vector must be three dimensional.');
                outputMatrix = [];
                return
            end
            % Apply element mapping
            outputMatrix = zeros(3,3);
            outputMatrix(1,2) = -inputVector(3);
            outputMatrix(1,3) =  inputVector(2);
            outputMatrix(2,1) =  inputVector(3);
            outputMatrix(2,3) = -inputVector(1);
            outputMatrix(3,1) = -inputVector(2);
            outputMatrix(3,2) =  inputVector(1); % Arrange the c
        end
        % CALL A UNIT TRIAD TO PLOT
        function [triadVectors] = drawTriad(position,R)
           colourVector = {'r','g','b'};
           triadVectors = eye(3); 
           for axis = 1:size(triadVectors,2)
               triadVectors(:,axis) = R*triadVectors(:,axis);
               % Draw vectors
%                quiver3(gca,position(1),position(2),position(3),triadVectors(1,axis),triadVectors(2,axis),triadVectors(3,axis),colourVector{axis});
           end
        end
        
        %% EULAR ROTATION TOOLS
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
        function [axisRates] = eularRatesToBodyAxisRates(eulerPosition,eulerRates)
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
           axisRates = conversionMatrix*eulerRates;
        end
        % EULER ROTATION MATICES (STANDARD AEROSPACE 'ZYX' ROTATION ORDER)
        function [R_BG,R_GB] = eulerToRotationMatrix(eulers)
            % This function converts the euler rotations (phi, theta, psi)
            % in the equivalent rotation matrix (R_AB) and inverse (R_BA).
            % "This matrix gives the conversion from body to navigation
            % frame"
            R_BG = [ cos(eulers(2))*cos(eulers(3)),...
                     cos(eulers(2))*sin(eulers(3)),...
                     -sin(eulers(2));...R
                     cos(eulers(3))*sin(eulers(1))*sin(eulers(2)) - cos(eulers(1))*sin(eulers(3)),...
                     cos(eulers(1))*cos(eulers(3)) + sin(eulers(1))*sin(eulers(3))*sin(eulers(2)),...
                     cos(eulers(2))*sin(eulers(1));...
                     sin(eulers(1))*sin(eulers(3)) + cos(eulers(1))*cos(eulers(3))*sin(eulers(2)),...
                     cos(eulers(1))*sin(eulers(3))*sin(eulers(2)) - cos(eulers(3))*sin(eulers(1)),...
                     cos(eulers(1))*cos(eulers(2))];                                % Correct for axis convention
            R_GB = transpose(R_BG);     
        end
        % GET THE ROTATION ANGLES BETWEEN TWO VECTORS
        function [angles,R]  = getVectorRotations(referenceVector,inputVector)
            % This function is designed to calculate the angles between a
            % given reference vector and an input vector; and a rotation
            % matrix between them.
            rotationAxis = cross(referenceVector,inputVector);
            rotationSkewMatrix = OMAS_axisTools.skew(rotationAxis);
            s = abs(sqrt(sum(rotationAxis.^2)));         % Sin of angle
            c = dot(referenceVector,inputVector);        % Cos of angle
            R = eye(3) + rotationSkewMatrix + rotationSkewMatrix^2*((1-c)/s^2);
            % Get angles
            angles = abs([0;asin(s);acos(c)]);
        end
        
        %% GENERIC QUATERNION TOOLS
        % ANALYTICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
        function [q] = getAnalyticalTriadRotation(referenceTriad,targetTriad)
            % This function defines the quaternion describing the rotations
            % between a reference triad and a second triad.
            
            % NORMALISE THE TRIAD VECTORS
            for dim = 1:size(referenceTriad,2)
                [referenceTriad(:,dim)] = referenceTriad(:,dim)/norm(referenceTriad(:,dim));
                [targetTriad(:,dim)] = targetTriad(:,dim)/norm(targetTriad(:,dim));
            end
            
            % EXTRACT REFERENCE AXES
            xAxis = targetTriad(:,1);
            zAxis = targetTriad(:,3);   % Get the rotated body axes
            xAxis_ref = referenceTriad(:,1);
            zAxis_ref = referenceTriad(:,3); % Get the reference unit ENU triad (trying to get to)
            
            % CORRECTION PLACEHOLDER
            R_correction = eye(3);
            
            % FIRST ALIGN THE Z-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GET THE QUATERNION ROTATION TO ALLIGN THE Z-AXES
            [q_zAlign] = OMAS_axisTools.vectorsToQuaternion(zAxis_ref,zAxis);    % Quaternion aligning global and body z vectors
            [R_zAlign] = OMAS_axisTools.quaternionToRotationMatrix(q_zAlign);    % Equivalent rotation matrix
            % ALIGN THE X AXIS IN THE Z-PLANE
            xAxis_intermediate = R_zAlign*xAxis;
            % TAKE ITS PROJECTIONS IN THE XY PLANE & RENORMALISE
            xAxis_intermediate(3) = 0;
            [xAxis_intermediate] = xAxis_intermediate/norm(xAxis_intermediate);
            % OTHERWISE JUST ALIGN THE X-AXIS VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GET THE QUATERNION ROTATION TO ALLIGN THE X-AXES
            [q_xAlign] = OMAS_axisTools.vectorsToQuaternion(xAxis_ref,xAxis_intermediate);
            [R_xAlign,~] = OMAS_axisTools.quaternionToRotationMatrix(q_xAlign);
            
            % COMPUTE THE COMPOSITE ROTATION MATRIX
            comp_rotation = R_xAlign * R_zAlign * R_correction;
            % COVERT THE ROTATION MATRIX TO A QUATERNION DESCRIBING THE
            % ROTATION FROM TRIAD_REF TO TRIAD_FINAL
            q = rotm2quat(comp_rotation);
        end
        % NUMERICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
%         function [q] = getNumericalTriadRotation(referenceTriad,targetTriad)
%            % This function defines the quaternion nescessary to rotate one 
%            % axis triad into alignment with a second triad. This is done by
%            % solving the rotation matrix-defined symbolic expressions for
%            % the equivalent quaternion rotation.
%            
%            syms q1 q2 q3 q4 real
%            
%            % Define the symbolic quaternion rotation matrix
%            R_q = zeros(3);
%            R_q(1,1) = q1^2 + q2^2 - q3^2 - q4^2;
%            R_q(1,2) = 2*(q1*q4 + q2*q3);
%            R_q(1,3) = 2*(q2*q4 - q1*q3);
%            R_q(2,1) = 2*(q3*q2 - q1*q4);
%            R_q(2,2) = q1^2 - q2^2 + q3^2 - q4^2;
%            R_q(2,3) = 2*(q1*q2 + q3*q4);
%            R_q(3,1) = 2*(q1*q3 + q2*q4);
%            R_q(3,2) = 2*(q3*q4 - q1*q2);
%            R_q(3,3) = q1^2 - q2^2 - q3^2 + q4^2;
%            % DEFINE THE SYMBOLIC EXPRESSIONS
%            expressions = [R_q*referenceTriad(:,1) - targetTriad(:,1);
%                           R_q*referenceTriad(:,2) - targetTriad(:,2);
%                           R_q*referenceTriad(:,3) - targetTriad(:,3);
%                           q1^2 + q2^2 + q3^2 + q4^2 - 1];
%            expressions = symfun(expressions,[q1 q2 q3 q4]);                % Define symbolic function in q_{1-4}
%            modeqn = @(x)double(expressions(x(1),x(2),x(3),x(4)));          % Rearrange to be in terms of q1,q2 etc
%            % DEFINE SOLVER OPTIONS
%            OPTIONS = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%            % PASS SYMBOLIC EQUATIONS TO F-SOLVE
%            x0 = [ 1 0 0 0];                 % INITIAL QUATERNION GUESS
%            q = fsolve(modeqn,x0,OPTIONS);   % SOLVE FOR THE ROTATION QUATERNION
%            % RE-NORMALISE THE QUATERNION
%            [q] = OMAS_axisTools.qUnit(q);
%         end
        % UPDATE A GIVEN QUATERNION
        function [q_new] = updateQuaternion(q_old,axisRates,dt)
            % This function is designed to compute the quaternion update routine from
            % one quaternion angle to all the appropriate representations.
            % INPUTS
            % q         - The initial quaternion axis reference
            % axisRates - The current axis rates [p;q;r]
            % dt        - The timestep between instances 
            % OUTPUTS
            % q         - The updated quaternion
            % q_dot     - The quaternion difference
            % R_q       - The rotation matrix R_G to R_B
            
            % Calculate the integration drift offset- lambda
            K = dt/dt; % k*dt <= 1          
            lambda = 1 - (q_old(1)^2 + q_old(2)^2 + q_old(3)^2 + q_old(4)^2);
            % Calculate the quaternion difference
            q_dot = 0.5*[K*lambda,  axisRates(3),  -axisRates(2),  axisRates(1);
                      -axisRates(3),    K*lambda,   axisRates(1),  axisRates(2);
                       axisRates(2), -axisRates(1),     K*lambda,  axisRates(3);
                      -axisRates(1), -axisRates(2),-axisRates(3),   K*lambda]*q_old;
            % Integrate the quaternion difference
            q_new = q_dot*dt + q_old;
            % Re-normalise
            q_new = OMAS_axisTools.qUnit(q_new);
        end
        % THE QUATERNION DIFFERENCE
        function [dq] = qDifference(q1,q2)
            % This function gets the quaternion that will rotate from the
            % q1 orientation to the q2 orientation [q1->q2].
            
            assert(size(q1,1) == 4,'The initial quaternion (q1) must be a 4x1 column vector')
            assert(size(q1,1) == 4,'The terminal quaternion (q2) must be a 4x1 column vector')
            
            % Get the quaternion inverse
            [q1_inv] = OMAS_axisTools.qInverse(q1);
            % Calculate difference
            dq = OMAS_axisTools.qMultiply(q1_inv,q2);
            % Re-normalise
            dq = OMAS_axisTools.qUnit(dq);
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
            q = OMAS_axisTools.qUnit(q); 
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
            [q] = OMAS_axisTools.qUnit(q);
        end
        % CONVERT ROTATION ANGLES INTO A QUATERNION
        function [q] = eulerToQuaternion(eulerAngles)
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
        function [eulerAngles] = quaternionToEuler(q)
            % Gets the rotation angles from an equivalent quaternion
            % attitude. Associated block:
            % "Quaternions to Rotation Angles"
            if length(q) ~= 4
                error('Input quaternion is of length %s',num2str(length(q)));
            end
            % Normalise the quaternion
            q = OMAS_axisTools.qUnit(q);
            eulerAngles = zeros(3,1);
            % Define the euler angles
            eulerAngles(1) = atan2((2*(q(2)*q(3) + q(1)*q(4))),(q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2));  % Phi   (roll)
            eulerAngles(2) = asin(-2*(q(2)*q(4) - q(1)*q(3)));                                        % Theta (pitch)
            eulerAngles(3) = atan2((2*(q(3)*q(4) + q(1)*q(2))),(q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2));  % Psi   (yaw)
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
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R_AB,R_BA] = quaternionToRotationMatrix(q)
            % This function defines the equivalent rotation matrix from the
            % provided quaternion. (Associated block "Quaternion to DCM")
            % INPUT:
            % q    - The quaternion rotation
            % OUTPUT:
            % R_AB - The rotation matrix through A to B
            % R_BA - The reverse rotation matrix from B to A
            
            assert(numel(q) == 4,'Input quaternion has wrong dimensions');
            
            R_AB = zeros(3,3);           
            % Normalise the quaternion
            q = OMAS_axisTools.qUnit(q);
            % Assemble the quaternion elements
            R_AB(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
            R_AB(1,2) = 2*(q(1)*q(4) + q(2)*q(3));
            R_AB(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
            R_AB(2,1) = 2*(q(3)*q(2) - q(1)*q(4));
            R_AB(2,2) = q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2;
            R_AB(2,3) = 2*(q(1)*q(2) + q(3)*q(4));
            R_AB(3,1) = 2*(q(1)*q(3) + q(2)*q(4));
            R_AB(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
            R_AB(3,3) = q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2;
            % Transpose the matrix for the global-body transformations
            R_BA = transpose(R_AB);
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
            qDiv = qDiv./(OMAS_axisTools.qNorm(r));
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
           qv = OMAS_axisTools.qUnit(qv);
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
            qUnit = q./OMAS_axisTools.qNorm(q);
        end
        % THE QUATERNION NORM
        function [qNorm] = qNorm(q)
            % GETS THE NORM OF THE QUATERNION VECTOR
            qNorm = sqrt(sum(q.^2));    % VALIDATED
        end
    end
end

