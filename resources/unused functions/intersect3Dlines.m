% INTERSECT 3D CARTESIAN VECTORS
function [P_intersect] = intersect3Dlines(PA,PB)
% Find intersection point of lines in 3D space, in the least squares sense.
% PA :          Nx3-matrix containing starting point of N lines
% PB :          Nx3-matrix containing end point of N lines
% P_Intersect : Best intersection point of the N lines, in least squares sense.
% distances   : Distances from intersection point to the input lines
% Anders Eikenes, 2012

Si = PB - PA; %N lines described as vectors
ni = Si ./ (sqrt(sum(Si.^2,2))*ones(1,3)); %Normalize vectors
nx = ni(:,1); ny = ni(:,2); nz = ni(:,3);
SXX = sum(nx.^2-1);
SYY = sum(ny.^2-1);
SZZ = sum(nz.^2-1);
SXY = sum(nx.*ny);
SXZ = sum(nx.*nz);
SYZ = sum(ny.*nz);
S = [SXX SXY SXZ;SXY SYY SYZ;SXZ SYZ SZZ];
CX  = sum(PA(:,1).*(nx.^2-1) + PA(:,2).*(nx.*ny)  + PA(:,3).*(nx.*nz));
CY  = sum(PA(:,1).*(nx.*ny)  + PA(:,2).*(ny.^2-1) + PA(:,3).*(ny.*nz));
CZ  = sum(PA(:,1).*(nx.*nz)  + PA(:,2).*(ny.*nz)  + PA(:,3).*(nz.^2-1));
C   = [CX;CY;CZ];
P_intersect = (S\C)';
end