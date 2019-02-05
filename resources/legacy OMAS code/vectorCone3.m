
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

% Iterate through the node
for nodeA = 1:nodes
    % GET ROTATED VECTOR
    [unit_hypVector] = rotateVectorAboutAxis(unit_axisVector,unit_hypVector,theta);
    % PLOT THE VECTORS
    q = quiver3(pointA(1),pointA(2),pointA(3),unit_hypVector(1),unit_hypVector(2),unit_hypVector(3),'g'); 
    q.AutoScaleFactor = 1;
    
    
    
    
    
    % GENERATE THE MESH
    for nodeB = 1:nodes
        % Get the point coordinate
        coordinate = (nodeB*nodalLength)*unit_hypVector;
        
        q = quiver3(pointA(1),pointA(2),pointA(3),coordinate(1),coordinate(2),coordinate(3),'g'); 
        q.AutoScaleFactor = 1;
        
        coordinateSet = horzcat(coordinateSet,coordinate);
%         xVector = horzcat(xVector,coordinate(1));
%         yVector = horzcat(yVector,coordinate(2));
%         zVector = horzcat(zVector,coordinate(3));
    end
end


cone = coordinateSet;




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