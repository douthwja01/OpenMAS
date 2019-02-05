% clear all; close all;

% INPUT VECTORS
% n = 20;
% colour = 'b';
% centerLineVector   = [0;10;0];    % An example rotation axis
% radiusVector = [0;-2;5];           % Radius of the cone

function [cone] = vectorCone2(centerLineVector,radiusVector,n) 
% This function is designed to construct a cone from two defining vectors.
% INPUTS:
% axisVector   - The centerline of the cone
% radiusVector - The radius of the cone
% 
% % PLOT CONFIGURATION
f = figure(1);
hold on;
grid on;
axis equal;
xlabel('x_{m}');
ylabel('y_{m}');
zlabel('z_{m}');

% % PLOT INITIAL PROBLEM
% q = quiver3(0,0,0,centerLineVector(1),centerLineVector(2),centerLineVector(3),'b'); 
% q.AutoScaleFactor = 1; 
% q = quiver3(centerLineVector(1),centerLineVector(2),centerLineVector(3),radiusVector(1),radiusVector(2),radiusVector(3),'r'); 
% q.AutoScaleFactor = 1; 
% 
% % GET HYPOTENUSE DATA
hypVector = centerLineVector + radiusVector;

% GET THE ROTATION ANGLES
unit_centerLineVector = centerLineVector/sqrt(sum(centerLineVector.^2));

theta = (2*pi/n);

for node = 1:n
    
    angle = node*theta;

    % GET ROTATED VECTOR
    [newVector] = rotateVectorAboutAxis(unit_centerLineVector,hypVector,angle);


    q = quiver3(0,0,0,newVector(1),newVector(2),newVector(3),'g'); 
    q.AutoScaleFactor = 1; 

end

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