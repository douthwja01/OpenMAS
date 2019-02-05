
function [cone] = vectorCone3(pointA,pointB,radialPoint,nodes) 
% This function is designed to generate a 3D cone mesh between two points
% and a radial definition using a certain number of nodes
% INPUTS:
% pointA - First 3D point
% pointB - Second 3D point
% radiuspoint - Radial point (cone desciption point) from point B


%% PLOT INITIAL PROBLEM
plot3(gca,pointA(1),pointA(2),pointA(3),'or'); 
plot3(gca,pointB(1),pointB(2),pointB(3),'og');
plot3(gca,radialPoint(1),radialPoint(2),radialPoint(3),'ob'); 

% CALCULATE THE VECTOR DEFINITIONS
axisVector = pointB - pointA;                          % Get the axisVector
unit_axisVector = axisVector/sqrt(sum(axisVector.^2)); % Unit axis vector
radialVector = radialPoint - pointB;                   % Calculate the radius vector
hypVector =  radialPoint - pointA;                     % Get the hypotenuseVector
mod_hypVector  = sqrt(sum(hypVector.^2));              % Modulus hyoptenuse
unit_hypVector = hypVector/mod_hypVector;              % Unit hypotenuse vector 

% GET VECTOR ROTATIONS
[angles,~] = getVectorRotations([1;0;0],axisVector);

% PLOT THE PROBLEM DEFINITIONS
q = quiver3(gca,pointB(1),pointB(2),pointB(3),radialVector(1),radialVector(2),radialVector(3),'g');
q.AutoScaleFactor = 1;
q = quiver3(gca,pointA(1),pointA(2),pointA(3),axisVector(1),axisVector(2),axisVector(3),'r'); 
q.AutoScaleFactor = 1;
q = quiver3(gca,pointA(1),pointA(2),pointA(3),hypVector(1),hypVector(2),hypVector(3),'b'); 
q.AutoScaleFactor = 1;

% CALCULATE THE CONE PROPERTIES
theta = (2*pi/nodes);               % Unit rotation angle
nodalLength = mod_hypVector/nodes;  % Unit length

% Containers
xVector = [];
yVector = [];
zVector = [];
coordinateSet = [];

%% BUILD THE CONE SURFACE
set1 = linspace(0,nodes,11);    % Build the X grid set (radially)
set2 = linspace(0,2*pi,nodes);  % Build the Y grid set (circular set)

set1 = 


[R,A] = meshgrid(set1,set2)     % Generate mesh grid 

m = nodalLength;
X = R .* cos(A) % Multiply Xgrid set by cos(Ygrid)
Y = R .* sin(A) % Multiply Ygrid set by sin(Ygrid)
Z = m*R; % Scale the z axis set

offset = pointA;
X = X + offset(1)
Y = Y + offset(2)
Z = Z + offset(3)

theta = real(angles(2)) 
phi = real(angles(3))

%[R,RX,RY,RZ] = getRotationMatrix(phi,theta,psi);

% Rotate in the azimuth
X1 = X*cos(phi) - Z*sin(phi);
Y1 = Y;
Z1 = X*sin(phi) + Z*cos(phi);
% Rotate in the elevation
X2 = X1*cos(theta) - Y1*sin(theta);
Y2 = X1*sin(theta) + Y1*cos(theta);
Z2 = Z1;

cone = mesh(X2,Y2,Z2);

% Iterate through the node
% for nodeA = 1:nodes
    % GET ROTATED VECTOR
%     [unit_hypVector] = rotateVectorAboutAxis(unit_axisVector,unit_hypVector,theta);
    % PLOT THE VECTORS
%     q = quiver3(pointA(1),pointA(2),pointA(3),unit_hypVector(1),unit_hypVector(2),unit_hypVector(3),'g'); 
%     q.AutoScaleFactor = 1;
    
    % GENERATE THE MESH
%     for nodeB = 1:nodes
%         % Get the point coordinate
%         coordinate = (nodeB*nodalLength)*unit_hypVector;
%         
%         q = quiver3(pointA(1),pointA(2),pointA(3),coordinate(1),coordinate(2),coordinate(3),'g'); 
%         q.AutoScaleFactor = 1;
%         
%         coordinateSet = horzcat(coordinateSet,coordinate);
% %         xVector = horzcat(xVector,coordinate(1));
% %         yVector = horzcat(yVector,coordinate(2));
% %         zVector = horzcat(zVector,coordinate(3));
%     end
% end


% cone = coordinateSet;

% cone = 0;

%cone = [xVector;yVector;zVector];
% URL=http://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit



% 
% % GET THE ROTATION ANGLES
% unit_centerLineVector = pointA/sqrt(sum(pointA.^2));
% 
% theta = (2*pi/n);
% 
% for node = 1:n
%     
%     angle = node*theta;
% 
%     % GET ROTATED VECTOR
%     [newVector] = rotateVectorAboutAxis(unit_centerLineVector,hypVector,angle);
% 
% 
%     q = quiver3(0,0,0,newVector(1),newVector(2),newVector(3),'g'); 
%     q.AutoScaleFactor = 1; 
% 
% end
% cone = 0
end
% VECTOR ROTATION FUNCTION
function [R,RX,RY,RZ] = getRotationMatrix(phi,theta,psi)
% This function is designed to rotate a given vector through the
% angles given in the angle vector.

% DECLARE THE INPUT ANGLES
%            phi   = angleVector(1); % rotation around the vector
%            theta = angleVector(2); % rotation in the xz plane (elevation)
%            psi   = angleVector(3); % rotation in the xy plane (azimuth)

% DEFINE THE EULAR ROTATON MATRIX
RX(1,1) = 1;
RX(1,2) = 0;
RX(1,3) = 0;
RX(2,1) = 0;
RX(2,2) = cos(phi);
RX(2,3) = -sin(phi);
RX(3,1) = 0;
RX(3,2) = sin(phi);
RX(3,3) = cos(phi);

RY(1,1) = cos(theta);
RY(1,2) = 0;
RY(1,3) = sin(theta);
RY(2,1) = 0;
RY(2,2) = 1;
RY(2,3) = 0;
RY(3,1) = -sin(theta);
RY(3,2) = 0;
RY(3,3) = cos(theta);

RZ(1,1) = cos(psi);
RZ(1,2) = -sin(psi);
RZ(1,3) = 0;
RZ(2,1) = sin(psi);
RZ(2,2) = cos(psi);
RZ(2,3) = 0;
RZ(3,1) = 0;
RZ(3,2) = 0;
RZ(3,3) = 1;

R = RZ*RY*RX;

% DEFINE THE INVERSE ROTATION MATRIX
Rinv = transpose(R);   % Invert for body to global
end
% GET THE ROTATIONS IN THE GLOBAL SYSTEM
function [angles,R] = getVectorRotations(referenceVector,inputVector)
% This function is designed to calculate the angles between a
% given reference vecto and an input vector and a rotation
% matrix between them.

rotationAxis = cross(referenceVector,inputVector);
rotationSkewMatrix = getSkewSymmetricMatrix(rotationAxis);
s = abs(sqrt(sum(rotationAxis.^2)));         % Sin of angle
c = dot(referenceVector,inputVector);        % Cos of angle
R = eye(3) + rotationSkewMatrix + rotationSkewMatrix^2*((1-c)/s^2);
% Get angles
angles = [0;asin(s);acos(c)];

end
% ROTATE VECTOR THROUGH AN ANGLE, AROUND A GIVEN AXIS
function [newVector] = rotateVectorAboutAxis(axisVector,oldVector,theta)
% This function is designed to calculate a vector
% following a rotation around a given axis vector, through a
% given angle.

% NORMALISE THE AXIS VECTOR
axisVector = axisVector/(sqrt(sum(axisVector.^2)));  % Normalize rotation axis

% GET THE CROSS PRODUCT PROECTION BETWEEN THE AXIS AND VECTOR
crossVector = cross(axisVector,oldVector);

% DETERMINE THE MAPPING OF EACH VECTOR COMPONENTS
newVector = cos(theta)*oldVector ...
    + (crossVector)*sin(theta)  ...
    + axisVector*(dot(axisVector,oldVector))*(1 - cos(theta));
end
% GENERATE SKEW-SYMMETRIC MATRIX
function [outputMatrix] = getSkewSymmetricMatrix(inputVector)
% This function generates a skew-symmetric cross-product matrix
% of the input vector. This allows a vector to be multiplied by
% the matrix and have the cross product be computed.
outputMatrix = zeros(3,3);
outputMatrix(1,2) = -inputVector(3);
outputMatrix(1,3) =  inputVector(2);
outputMatrix(2,1) =  inputVector(3);
outputMatrix(2,3) = -inputVector(1);
outputMatrix(3,1) = -inputVector(2);
outputMatrix(3,2) =  inputVector(1); % Arrange the components
end