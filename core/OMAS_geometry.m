%% OPENMAS 3D-AXIS MANIPULATION TOOLS (OMAS_geometry.m) %%%%%%%%%%%%%%%%%%%
% This class provides the axis manipulation utilities used by OpenMAS to
% evaluate object progression over time.

% Author: James A. Douthwaite 06/09/2017

classdef OMAS_geometry
    % This class is designed to handle all the axis transformations and
    % conversions as part of the OpenMAS simulator.
    
    %% GENERIC TOOLS
    methods (Static) 
        % SIGNED ANGLE BETWEEN TO VECTORS
        function [planarAngle,isRightLogical] = signedPlanarAngle(n,u,v)
            % This function computes the signed angle between two vectors
            % with respect to a provided normal.
            [planarAngle,isRightLogical] = GetSignedPlanarAngle(n,u,v);
        end
        % RODRIGUES ROTATION
        function [v_rotated] = rodriguesRotation(u,k,theta)
            v_rotated = rodriguesRotation(u,k,theta);
        end
        % GET VECTOR PLANAR PROJECTION
        function [v_proj] = vectorPlanarProjection(n,v)
            v_proj = GetPlanarProjection(n,v);
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
        end
        % SKEW SYMMETRIC MATRIX
        function [vSkew] = skew(v)
            vSkew = skew(v);
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
            % Call the assoicated mex file
            [eulerRates] = GetEularRatesFromOmega(eulerPosition,bodyAxisRates);
        end
        % CONVERT EULER RATES TO BODY AXIS RATES
        function [bodyAxisRates] = eularRatesToBodyAxisRates(eulerPosition,eulerRates)
            % This function calculates the conversion between the eular
            % rates and the body axis rates.
            % INPUT HANDLING
            assert(numel(eulerRates) == 3,'Incorrect number of Euler rates');
            assert(numel(eulerPosition) == 3,'Incorrect number of Euler rotations');
            % Call associated mex file
            bodyAxisRates = GetOmegaFromEulerRates(eulerPosition,eulerRates);
        end
        % EULER ROTATION MATRIX TO EULER ANGLES
        function [eulers] = rotationMatrixToEulers(R)
            % Rotation of the fixed to the global pose
            eulers = GetEulersFromRotationMatrix(R);
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
            % Get the mex rotations
            R = R_eta(eulers);
        end
        % EULER ROTATION MATRIX (Roll about the X-axis)
        function [R_phi]   = eulerToRotationMatrix_roll(phi)
            % This rotation matrix maps a vector through the rotation angle
            % 'phi' about a reference 'x' axis.
            % ROLL AXIS ROTATIONS
            R_phi = R_x(phi);
        end
        % EULER ROTATION MATRIX (Pitch about the Y-axis)
        function [R_theta] = eulerToRotationMatrix_pitch(theta)
            % This rotation matrix maps a vector through the rotation angle
            % 'theta' about a reference 'y' axis.
            % PITCH AXIS ROTATIONS
            R_theta = R_y(theta);
        end
        % EULER ROTATION MATRIX (Yaw about the Z-axis)
        function [R_psi]   = eulerToRotationMatrix_yaw(psi)
            % This rotation matrix maps a vector through the rotation angle
            % 'psi' about a reference 'z' axis.
            % YAW AXIS ROTATIONS
            R_psi = R_z(psi);
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
%             q = GetAnalyticalTriadRotation(referenceTriad,targetTriad);
            q = GetAnalyticalTriadRotation(referenceTriad,targetTriad);
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
            [q_dot] = qDifferential(q0,omega);    
            % Integrate the quaternion difference
            q_new = q_dot*dt + q0;
            % Re-normalise
            q_new = unit(q_new);
        end            
        % GET THE QUATERNION DIFFERNTIAL (as a function of omega)
        function [q_dot] = quaternionDifferential(q0,omega)
            % This function was created to allow the quaternion
            % differential to be called seperately to in the integration
            % function. This allows the quaternion 'q0' to be integrated
            % externally.
            q_dot = qDifferential(q0,omega);
        end
        % ROTATE A VECTOR THROUGH A QUATERNION
        function [v1] = quaternionVectorRotation(q,v0)
            % This function rotates a cartesian vector through a quaternion
            % rotation. Associated block:
            % "Quaternion Rotation"
            v1 = qVectorRotation(q,v0);
        end
        % GET THE QUATERNION BETWEEN TWO VECTORS
        function [q] = vectorsToQuaternion(u,v)
            q = qArgument(u,v);
        end
        % CONVERT ROTATION ANGLES INTO A QUATERNION
        function [q] = eulersToQuaternion(eulerAngles)
            % THis function converts euler angle rotations into a
            % quaternion via its components. Associated block:
            % "Rotation Angles to Quaternions".
            q = GetQuaternionFromEulers(eulerAngles);
        end
        % CONVERT QUATERNION INTO ROTATION ANGLES
        function [eta] = quaternionToEulers(q)  % NOTATION B
            % Gets the rotation angles from an equivalent quaternion
            % attitude. Associated block:
            % "Quaternions to Rotation Angles"
            assert(size(q,1) == 4,'A quaternion column vector is expected 4x1');
            % Call the associated mex file
            eta = qRotations(q);
        end
        % CONVERT ROTATION MATRIX TO QUATERNION
        function [q] = rotationMatrixToQuaternion(R)
            % This function is designed to convert from a rotation matrix
            % to an equivalent quaternion. This function is also parallel
            % to "rotm2quat.m".
            % Call the associated mex file
            q = GetQuaternionFromRotationMatrix(R);
        end
        % GET ROTATION MATRIX FROM QUATERNION
        function [R] = quaternionToRotationMatrix(q)
            % This function defines the equivalent rotation matrix from the
            % provided quaternion. (Associated block "Quaternion to DCM")
            % INPUT:
            % q    - The quaternion rotation
            R = R_q(q);
%             R = quat2rotm(q');  % Call the robotics toolbox
        end
        % THE QUATERNION DIFFERENCE
        function [dq] = qDifference(q,v)
            % This function gets the quaternion that will rotate from the
            % q1 orientation to the q2 orientation [q1->q2].
            dq = qDifference(q,v);                                     % Call the associated mex file
        end       
        % QUATERNION DIVISION
        function [qDiv] = qDivide(q,v)
            % The function compute the division of quaternion q, by
            % quaterion v.
            % Associated block:
            % "Quaternion Division"
            qDiv = qDivide(q,v);                                       % Call the associated mex file
        end
        % QUATERNION MULTIPLICATION
        function [qv] = qMultiply(q,v)
            % Calculate the product of two quaternions
            % Associated block:
            % "Quaternion Multiplication"
            % Multiply the quaternion elements
            qv = qMultiply(q,v);                                       % Call the associated mex file
        end
        % QUATERNION INVERSE
        function [qInv] = qInverse(q)
            % The quaternion norm
            qInv = qInverse(q);                                        % Call the associated mex file
        end
        % QUATERNION CONJUGATE
        function [qCon] = qConjugate(q)
            % This function defines the conjugate of the given quaternion.
            % Associated block:
            % "Quaternion Conjugate"
            qCon = qConjugate(q);
        end
    end
    %% 3D GEOMETRIC TOOLS
    methods (Static)
        % INTERSECT - SPHERE AND ROTATED CUBOID (OBB)
        function [intersectFlag] = intersect_OBB_sphereCuboid(centerA,radiusA,centerB,cuboidVerticesB)
            % This function computes the intersection between a sphere and a 
            % rotated bounding box (OBB). This function assumes that the
            % problem is being resolved in the axes of A (i.e) 'R' is the
            % rotation matrix taking the cube from its own frame to the
            % axes of A.
            intersectFlag = CheckOBBSphereIntersection(centerA,radiusA,centerB,cuboidVerticesB);
        end
        % INTERSECT - TWO ROTATED CUBOIDS (OBB)
        function [intersectFlag] = intersect_OBB_cuboids(centerA,cuboidVerticesA,centerB,cuboidVerticesB)
            intersectFlag = CheckOBBOBBIntersection(centerA,cuboidVerticesA,centerB,cuboidVerticesB);
        end
        % INTERSECT - SPHERE AND AXIS ALIGNED CUBOID (AABB)
        function [intersectFlag] = intersect_AABB_sphereCuboid(center,radius,cubeMins,cubeMaxs)
            % get cuboid closest point to sphere center by clamping
            %             x = max(cubeMins(1), min(sphere.x, cubeMaxs(1)));
            %             y = max(cubeMins(2), min(sphere.y, cubeMaxs(2)));
            %             z = max(cubeMins(3), min(sphere.z, cubeMaxs(3)));
            intersectFlag = CheckAABBSphereIntersection(center,radius,cubeMins,cubeMaxs);
        end
        % INTERSECT - AXIS ALIGNED CUBOIDS (AABB)
        function [intersectFlag] = intersect_AABB_cuboids(minA,maxA,minB,maxB)
            % This function determines if two 3D cuboids currently intersect
            intersectFlag = CheckAABBAABBIntersection(minA,maxA,minB,maxB);
        end
        % INTERSECT - SPHERES
        function [intersectFlag] = intersect_spheres(positionA,radiusA,positionB,radiusB)
            % This uses a simple one dimensional comparison to determine
            % the overlap of two spherical volumes.
            intersectFlag = CheckSphereSphereIntersection(positionA,radiusA,positionB,radiusB);
        end
        % INTERSECT - RANGE
        function [intersectFlag] = rangeIntersect(minA,maxA,minB,maxB)
            % This function determines the intersection between two 1D
            % ranges, and the distance between them.
            intersectFlag = CheckRangeIntersection(minA,maxA,minB,maxB);
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
           
           % Call external mex function
           p_int = GetRaySphereIntersection(ray,centroid,radius);
        end
        % PROJECTION OF A RAY ON A PLANE DEFINED BY A NORMAL
        function [v_prj] = rayProjectionOnPlane(ray,planeNormal)
            % This function computes the projection of a vector on a plane
            % defined by its normal vector.
            v_prj = GetRayProjectionOnPlane(ray,planeNormal);
        end
        % PROJECTION OF A POINT ON A VECTOR
        function [v_mag] = pointProjectionOnRay(ray,point)
            % Calcuate the projection magnitude(distance)
            v_mag = GetPointProjectionOnRay(ray,point);
        end
        % DEFINE RAY
        function [rayObj] = defineRay(origin,direction,magnitude)
           % Define a ray with a given a origin and direction.
           assert(numel(origin) == 3 && numel(direction) == 3,'Parameters must be 3D');
           if nargin < 3
               magnitude = 1;
           end
           rayObj = ray(origin,direction,magnitude);
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

