% This function is designed to determine whether a sphere collides with a
% give triangle

function [haveCollided] = OMAS_sphereTriangleIntersection(P,r,A,B,C)

% INPUT CHECK
assert(numel(P) == 3,'Centroid must be a column vector [3x1]');
assert(numel(r) == 1,'Radius must be a scalar');
assert(numel(A) == 3 && numel(B) == 3 && numel(C) == 3,...
       'Vertex point must be a column vector [3x1]');

% CHECK PLANAR SEPARATION
A = A - P;
B = B - P;
C = C - P;
rr = r^2;
V = cross(B-A,C-A);
d = dot(A,V);
e = dot(V,V);
sep1 = d^2 > rr*e;
% CHECK VERTICES
aa = dot(A,A);
ab = dot(A,B);
ac = dot(A,C);
bb = dot(B,B);
bc = dot(B,C);
cc = dot(C,C);
sep2 = (aa > rr) & (ab > aa) & (ac > aa);
sep3 = (bb > rr) & (ab > bb) & (bc > bb);
sep4 = (cc > rr) & (ac > cc) & (bc > cc);
% CHECK TRIANGLE EDGES
AB = B - A;
BC = C - B;
CA = A - C;
d1 = ab - aa;
d2 = bc - bb;
d3 = ac - cc;
e1 = dot(AB,AB);
e2 = dot(BC,BC);
e3 = dot(CA,CA);
Q1 = A*e1 - d1*AB;
Q2 = B*e2 - d2*BC;
Q3 = C*e3 - d3*CA;
QC = C*e1 - Q1;
QA = A*e2 - Q2;
QB = B*e3 - Q3;
sep5 = dot(Q1,Q1) > (rr*e1*e1) & dot(Q1,QC) > 0;
sep6 = dot(Q2,Q2) > (rr*e2*e2) & dot(Q2,QA) > 0;
sep7 = dot(Q3,Q3) > (rr*e3*e3) & dot(Q3,QB) > 0;
isSeparated = sep1 || sep2 || sep3 || sep4 || sep5 || sep6 || sep7;
haveCollided = ~isSeparated;
end