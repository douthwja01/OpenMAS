
function pointCloud = defineCuboid(minsA,maxsA)
% Define the point cloud from the 

pointCloud = zeros(8,3);
pointCloud(1,:) = [maxsA(1),maxsA(2),maxsA(3)];
pointCloud(2,:) = [maxsA(1),maxsA(2),minsA(3)];
pointCloud(3,:) = [maxsA(1),minsA(2),minsA(3)];
pointCloud(4,:) = [maxsA(1),minsA(2),maxsA(3)];
pointCloud(5,:) = [minsA(1),minsA(2),minsA(3)];
pointCloud(6,:) = [minsA(1),minsA(2),maxsA(3)];
pointCloud(7,:) = [minsA(1),maxsA(2),maxsA(3)];
pointCloud(8,:) = [minsA(1),maxsA(2),minsA(3)];
end