
close all; clear all;


% interiorPoints = rand(200,3);      %# Generate 200 3-D points
% DT = DelaunayTri(interiorPoints);  %# Create the tetrahedral mesh
% hullFacets = convexHull(DT);       %# Find the facets of the convex hull
% 
% %# Plot the scattered points:
% subplot(2,2,1);
% scatter3(interiorPoints(:,1),interiorPoints(:,2),interiorPoints(:,3),'.');
% axis equal;
% title('Interior points');
% 
% %# Plot the tetrahedral mesh:
% subplot(2,2,2);
% tetramesh(DT);
% axis equal;
% title('Tetrahedral mesh');
% 
% %# Plot the 3-D convex hull:
% subplot(2,2,3);
% trisurf(hullFacets,DT.X(:,1),DT.X(:,2),DT.X(:,3),'FaceColor','c')
% axis equal;
% title('Convex hull');

% DT

% Create a set of points P
% [x1, y1, z1] = sphere(24);
% x1 = x1(:);
% y1 = y1(:);
% z1 = z1(:);
% x2 = x1+5;
% P = [x1 y1 z1; x2 y1 z1];
% P = unique(P,'rows');
% % Plot the points
% plot3(P(:,1),P(:,2),P(:,3),'.')
% axis equal
% % Use alphaShape to create a polyhedron that envelops the points
% figure
% shp = alphaShape(P,1)
% plot(shp,'FaceColor','cyan','FaceAlpha',0.3);


% Create a set of points (x,y)
th = (pi/12:pi/12:2*pi)';
x1 = [reshape(cos(th)*(1:5), numel(cos(th)*(1:5)),1); 0];
y1 = [reshape(sin(th)*(1:5), numel(sin(th)*(1:5)),1); 0];
x = [x1; x1+15; 25; 26; 26];
y = [y1; y1; 0; 0; 0.25];
shp = alphaShape(x,y,2)
plot(shp,'EdgeColor','none')
hold on
plot(x,y,'.')
hold off
% Query the areas of the individual regions
area(shp,1:numRegions(shp))
% Select a RegionThreshold of 0.2 to remove the small region
shp.RegionThreshold = 0.2
% Plot to observe the region removed
figure
plot(shp,'EdgeColor','none')