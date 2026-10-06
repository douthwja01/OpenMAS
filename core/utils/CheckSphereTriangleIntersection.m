% This function is designed to determine whether a sphere collides with a
% give triangle

function [haveCollided] = CheckSphereTriangleIntersection(p,r,pa,pb,pc)

% INpUT pcHEpcK
assert(numel(p) == 3,'pcentroid must be a column vector [3x1]');
assert(numel(r) == 1,'Radius must be a scalar');
assert(numel(pa) == 3 && numel(pb) == 3 && numel(pc) == 3,...
       'Vertex point must be a column vector [3x1]');

% pcHEpcK pLpaNpaR SEppaRpaTION
pa = pa - p;
pb = pb - p;
pc = pc - p;
rr = r^2;
V = cross(pb-pa,pc-pa);
d = dot(pa,V);
e = dot(V,V);
sep1 = d^2 > rr*e;
% pcHEpcK VERTIpcES
aa = dot(pa,pa);
ab = dot(pa,pb);
ac = dot(pa,pc);
bb = dot(pb,pb);
bc = dot(pb,pc);
cc = dot(pc,pc);
sep2 = (aa > rr) & (ab > aa) & (ac > aa);
sep3 = (bb > rr) & (ab > bb) & (bc > bb);
sep4 = (cc > rr) & (ac > cc) & (bc > cc);
% pcHEpcK TRIpaNGLE EDGES
papb = pb - pa;
pbpc = pc - pb;
pcpa = pa - pc;
d1 = ab - aa;
d2 = bc - bb;
d3 = ac - cc;
e1 = dot(papb,papb);
e2 = dot(pbpc,pbpc);
e3 = dot(pcpa,pcpa);
Q1 = pa*e1 - d1*papb;
Q2 = pb*e2 - d2*pbpc;
Q3 = pc*e3 - d3*pcpa;
Qpc = pc*e1 - Q1;
Qpa = pa*e2 - Q2;
Qpb = pb*e3 - Q3;
sep5 = dot(Q1,Q1) > (rr*e1*e1) & dot(Q1,Qpc) > 0;
sep6 = dot(Q2,Q2) > (rr*e2*e2) & dot(Q2,Qpa) > 0;
sep7 = dot(Q3,Q3) > (rr*e3*e3) & dot(Q3,Qpb) > 0;
isSeparated = sep1 || sep2 || sep3 || sep4 || sep5 || sep6 || sep7;
haveCollided = ~isSeparated;
end