
oldVector = [1;0;0];
axisVector = [0;0;1];
alpha = deg2rad(45);
display(rad2deg(alpha));

%theta = alpha - (pi/2);

[newVector] = rotateVectorAboutAxis(oldVector,axisVector,alpha)

figure(1);
hold on;
q = quiver3(0,0,0,axisVector(1),axisVector(2),axisVector(3),'b','filled');
q.AutoScaleFactor = 1;
q = quiver3(0,0,0,oldVector(1),oldVector(2),oldVector(3),'r','filled');
q.AutoScaleFactor = 1;
q = quiver3(0,0,0,newVector(1),newVector(2),newVector(3),'g','filled');
q.AutoScaleFactor = 1;

