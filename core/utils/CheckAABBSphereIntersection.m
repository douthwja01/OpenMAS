% INTERSECT - SPHERE AND AXIS ALIGNED CUBOID (AABB)
function [intersectFlag,separation_AABB] = CheckAABBSphereIntersection(center,radius,cubeMins,cubeMaxs)
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