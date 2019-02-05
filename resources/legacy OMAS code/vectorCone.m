% CONSTRUCT CONE FROM VECTOR DEFINTIONS
function [Cone] = vectorCone(pointA,pointB,radialPoint,nodes)
% This function is designed to construct a cone mesh between two 
% 3D points Using a radial Point to specify the radial properties.
% INPUTS:
% pointA - First cone axis point (origin).
% pointB - Second cone axis point.
% radialPoint - Third point in space that defines the cones maximum radius
% nodes  - The number of mesh nodes
% OUTPUTS
% Cone   - Cone mesh structure

% INPUT HANDLING
if ~exist('nodes','var')
    nodes = 10;
end
if ~exist('pointA','var')
    pointA = [0;0;0];
end

% MANUAL INPUTS
coneColour = 'g';
coneAlpha = 0.5;
coneClosed = 1;
coneLines = 1;

% DEFINE VECTOR PROBLEM PARAMETERS
axisVector = pointB - pointA;
mod_axisVector = sqrt(sum((axisVector).^2)); % Axis vector properties
tangent = radialPoint - pointA;
mod_tangent = sqrt(sum((tangent).^2));       % Tangental vector properties

% DEFINE THE TANGENT-AXIS PROJECTION
trueAB = (dot(tangent,axisVector)/mod_axisVector^2)*axisVector;
mod_trueAB = sqrt(sum((trueAB).^2));

% GET THE RADIUS MODULUS
mod_radius = sqrt(mod_tangent^2-mod_trueAB^2);

% Creating 2 circles in the YZ plane
t=linspace(0,2*pi,nodes)';      % Create axis point set
xa2 = zeros(length(t),1);
xa3 = zeros(size(xa2));
xb2 = mod_radius*cos(t);
xb3 = mod_radius*sin(t);        % Scale the second axis in the

% Creating the points in the X-Direction
x1=[0 mod_trueAB];

% Creating (Extruding) the cylinder points in the X-Directions
xx1 = repmat(x1,length(xa2),1);
xx2 = [xa2 xb2];
xx3 = [xa3 xb3]; % Concatinate second circle set

% Drawing two filled cirlces to close the cylinder
if coneClosed == 1
    EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
    EndPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'r');
end

% GENERATE THE CONE
% Plot the cone from the origin along x-axis and scale to size
Cone = mesh(xx1,xx2,xx3);

% Get the planar rotation angle
unit_Vx=[1 0 0];
angle_X1X2 = acos(dot(unit_Vx,axisVector)/(norm(unit_Vx)*mod_axisVector))*180/pi;

% Get rotation axis
axis_rot = cross([1 0 0],axisVector);

% Rotating the plotted Cone and the end plate circles to the required
% angles
if angle_X1X2~=0 % Rotation is not needed if required direction is along X
    rotate(Cone,axis_rot,angle_X1X2,[0 0 0])
    if coneClosed==1
        rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
        rotate(EndPlate2,axis_rot,angle_X1X2,[0 0 0])
    end
end

% Till now Cone has only been aligned with the required direction, but
% position starts from the origin. so it will now be shifted to the right
% position
if coneClosed == 1
    set(EndPlate1,'XData',get(EndPlate1,'XData') + pointA(1))
    set(EndPlate1,'YData',get(EndPlate1,'YData') + pointA(2))
    set(EndPlate1,'ZData',get(EndPlate1,'ZData') + pointA(3))
    
    set(EndPlate2,'XData',get(EndPlate2,'XData') + pointA(1))
    set(EndPlate2,'YData',get(EndPlate2,'YData') + pointA(2))
    set(EndPlate2,'ZData',get(EndPlate2,'ZData') + pointA(3))
end
set(Cone,'XData',get(Cone,'XData') + pointA(1))
set(Cone,'YData',get(Cone,'YData') + pointA(2))
set(Cone,'ZData',get(Cone,'ZData') + pointA(3))

% SET THE COLOUR OF THE CONE AND END PLATES
set(Cone,'AmbientStrength',1,...
    'FaceColor',coneColour,...
    'FaceLighting','gouraud',...
    'FaceAlpha',coneAlpha);        % Cone verticies
if coneClosed==1
    set([EndPlate1 EndPlate2],...
        'AmbientStrength',1,...
        'FaceColor',coneColour,...
        'FaceLighting','gouraud',...
        'FaceAlpha',coneAlpha);         % End-plate
else
    EndPlate1=[];
    EndPlate2=[];
end

% If lines are not needed making it disapear
if coneLines == 0
    set(Cone,'EdgeAlpha',0)
end
end