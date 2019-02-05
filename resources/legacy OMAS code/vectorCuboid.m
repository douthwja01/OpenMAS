
clear all; close all;

minVector = [-4;-1;-3];
maxVector = [2;4;5];
nodes = 20;
closed = 1;

% The unit vertices of the 3D feasible region
selectionMatrix = [0 1 0 0 1 1 0 1;...
                   0 0 1 0 0 1 1 1;...
                   0 0 0 1 1 0 1 1];
maxVector = repmat((maxVector-minVector),1,8);
FRmatrix = selectionMatrix.*maxVector + repmat(minVector,1,8); % Assemble points set

[k,v] = boundary(FRmatrix(1,:)',FRmatrix(2,:)',FRmatrix(3,:)')
trisurf(k,FRmatrix(1,:)',FRmatrix(2,:)',FRmatrix(3,:)','Facecolor','red','FaceAlpha',0.1)


T = delaunayn(k)
% tetramesh(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the point-relative geometry
% geometry = maxVector - minVector; 

% % Calculating the length of the Cone
% length_cyl = 4;
% xdim = geometry(1);
% ydim = geometry(2);
% zdim = geometry(3);
% 
% R = [3 4];
% 
% t = linspace(0,xdim,nodes)'
% xa2=R(1)*cos(t);
% xa3=R(1)*sin(t);
% xb2=R(2)*cos(t);
% xb3=R(2)*sin(t);
% 
% % Creating the points in the X-Direction
% x1=[0 length_cyl];
% 
% % Creating (Extruding) the cylinder points in the X-Directions
% xx1=repmat(x1,length(xa2),1);
% xx2=[xa2 xb2];%xx2=repmat(x2,1,2);
% xx3=[xa3 xb3];%xx3=repmat(x3,1,2);
% 
% % Drawing two filled cirlces to close the cylinder
% if closed==1
%     hold on
%     EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
%     EndPlate2=fill3(xx1(:,2),xx2(:,2),xx3(:,2),'r');
% end
% 
% % Plotting the cylinder along the X-Direction with required length starting
% % from Origin
% figure(1);
% Cone=mesh(xx1,xx2,xx3);
% 
% figure(2);
% grid on;

% Multiply-connected polygon - a square with a square hole.
% % Counterclockwise outer loop, clockwise inner loop.
% xv = [0 3 3 0 0 NaN 1 1 2 2 1]; %x vector coordinates
% yv = [0 0 3 3 0 NaN 1 2 2 1 1];
% 
% x = rand(1000,1)*3; y = rand(1000,1)*3; % Point cloud
% 
% 
% in = inpolygon(x,y,xv,yv);
% plot(xv,yv,x(in),y(in),'.r',x(~in),y(~in),'.b')




