%*******************************************************************************
% function:	p_poly_dist
% Description:	distance from point to polygon whose vertices are specified by the
%              vectors xv and yv
% Input:  
%    x - point's x coordinate
%    y - point's y coordinate
%    xv - vector of polygon vertices x coordinates
%    yv - vector of polygon vertices x coordinates
% Output: 
%    d - distance from point to polygon (defined as a minimal distance from 
%        point to any of polygon's ribs, positive if the point is outside the
%        polygon and negative otherwise)
%    x_poly: x coordinate of the point in the polygon closest to x,y
%    y_poly: y coordinate of the point in the polygon closest to x,y
%
% Routines: p_poly_dist.m
% Revision history:
%    03/31/2008 - return the point of the polygon closest to x,y
%               - added the test for the case where a polygon rib is 
%                 either horizontal or vertical. From Eric Schmitz.
%               - Changes by Alejandro Weinstein
%    7/9/2006  - case when all projections are outside of polygon ribs
%    23/5/2004 - created by Michael Yoshpe 
% Remarks:
%*******************************************************************************
function [d,x_poly,y_poly] = p_poly_dist(x, y, xv, yv) 

% If (xv,yv) is not closed, close it.
xv = xv(:);
yv = yv(:);
Nv = length(xv);
if ((xv(1) ~= xv(Nv)) || (yv(1) ~= yv(Nv)))
    xv = [xv ; xv(1)];
    yv = [yv ; yv(1)];
%     Nv = Nv + 1;
end

% linear parameters of segments that connect the vertices
% Ax + By + C = 0
A = -diff(yv);
B =  diff(xv);
C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);

% find the projection of point (x,y) on each rib
AB = 1./(A.^2 + B.^2);
vv = (A*x+B*y+C);
xp = x - (A.*AB).*vv;
yp = y - (B.*AB).*vv;

% Test for the case where a polygon rib is 
% either horizontal or vertical. From Eric Schmitz
id = find(diff(xv)==0);
xp(id)=xv(id);
clear id
id = find(diff(yv)==0);
yp(id)=yv(id);

% find all cases where projected point is inside the segment
idx_x = (((xp>=xv(1:end-1)) & (xp<=xv(2:end))) | ((xp>=xv(2:end)) & (xp<=xv(1:end-1))));
idx_y = (((yp>=yv(1:end-1)) & (yp<=yv(2:end))) | ((yp>=yv(2:end)) & (yp<=yv(1:end-1))));
idx = idx_x & idx_y;

% distance from point (x,y) to the vertices
dv = sqrt((xv(1:end-1)-x).^2 + (yv(1:end-1)-y).^2);

if(~any(idx)) % all projections are outside of polygon ribs
   [d,I] = min(dv);
   x_poly = xv(I);
   y_poly = yv(I);
else
   % distance from point (x,y) to the projection on ribs
   dp = sqrt((xp(idx)-x).^2 + (yp(idx)-y).^2);
   [min_dv,I1] = min(dv);
   [min_dp,I2] = min(dp);
   [d,I] = min([min_dv min_dp]);
   if I==1, %the closest point is one of the vertices
       x_poly = xv(I1);
       y_poly = yv(I1);
   elseif I==2, %the closest point is one of the projections
       idxs = find(idx);
       x_poly = xp(idxs(I2));
       y_poly = yp(idxs(I2));
   end
end

if(inpolygon(x, y, xv, yv)) 
   d = -d;
end