function [velocityA,Vel_Collide_union] = Avoid_Algorithm_3D_allRotate(positionA,positionB,velocityA,velocityB,Size,prev_collisionVelocity,currentVelocity,colour,precision)
%% Avoid Algorithm

%Positions of Points Shifted by VelB
Posa = [positionA(1)+velocityB(1),positionA(2)+velocityB(2),positionA(3)+velocityB(3),1];
Posb = [positionB(1)+velocityB(1),positionB(2)+velocityB(2),positionB(3)+velocityB(3),1];



Vel_No_Collide = ones(30000,4)*10;                        
Vel_Collide = ones(2000,4)*10;

Lambda_AB = sqrt((Posb(1) - Posa(1))^2 + (Posb(2) - Posa(2))^2 + (Posb(3) - Posa(3))^2);
if Lambda_AB <= 10
    Quad1_yz = 0;
    Quad2_yz = 0;
    Quad3_yz = 0;
    Quad4_yz = 0;
    x_b = Posb(1);
    %Location of Agent B xz plane
    if (positionB(1) > positionA(1)) && (positionB(3) > positionA(3))
        x_b = x_b + 3;
    elseif (positionB(1) > positionA(1)) && (positionB(3) <= positionA(3))
        x_b = x_b + 3;
    elseif (positionB(1) < positionA(1)) && (positionB(3) < positionA(3))
        x_b = x_b - 3;
    elseif (positionB(1) < positionA(1)) && (positionB(3) >= positionA(3))
        x_b = x_b - 3;
    end
    %Location of Agent B yz plane
    if (positionB(2) > positionA(2)) && (positionB(3) > positionA(3))   
        Quad1_yz = 1;
    elseif (positionB(2) > positionA(2)) && (positionB(3) <= positionA(3))
        Quad2_yz = 1;
    elseif (positionB(2) < positionA(2)) && (positionB(3) < positionA(3))
        Quad3_yz = 1;
    elseif (positionB(2) < positionA(2)) && (positionB(3) >= positionA(3))
        Quad4_yz = 1;
    else %(PosB(2) == PosA(2)) && (PosB(3) == PosA(3))
        Posb = Posb + 0.0000000000000001; 
    end

    if positionB(1) == positionA(1) && positionB(2) == positionA(2)
        if positionB(3) > positionA(3)
            z_b = positionB(3) + 3;
        elseif positionB(3) < positionA(3)
            z_b = positionB(3) - 3;
        end
        x_b = positionB(1);
        y_b = positionB(2);
        Posb = [x_b,y_b,z_b,1];
    elseif positionB(1) == positionA(1)
        if Quad1_yz == 1
            x_b = Posb(1);
            y_b = positionB(2) + 3;
        elseif Quad2_yz == 1
            x_b = positionB(1);
            y_b = positionB(2) + 3;
        elseif Quad3_yz == 1
            x_b = positionB(1);
            y_b = positionB(2) - 3;
        elseif Quad4_yz == 1
            x_b = positionB(1);
            y_b = positionB(2) - 3;
        end
        Grad_AB_y = (Posa(3) - Posb(3)) / (Posa(2) - Posb(2));
        c_yz = Posa(3) - Grad_AB_y*Posa(2);
        z_b = Grad_AB_y*y_b + c_yz;
        Posb = [x_b,y_b,z_b,1];
    elseif positionB(2) == positionA(2)
        y_b = Posb(2);
        Grad_AB_z = (Posa(3) - Posb(3)) / (Posa(1) - Posb(1));
        c_xz = Posa(3) - Grad_AB_z*Posa(1);
        z_b = Grad_AB_z*x_b + c_xz;
        Posb = [x_b,y_b,z_b,1];
    elseif positionB(3) == positionA(3)
        z_b = Posb(3);
        Grad_AB_xy = (Posa(2) - Posb(2)) / (Posa(1) - Posb(1));
        c_xy = Posa(2) - Grad_AB_xy*Posa(1);
        y_b = Grad_AB_xy*x_b + c_xy;
        Posb = [x_b,y_b,z_b,1];
    else
        Grad_AB_z = (Posa(3) - Posb(3)) / (Posa(1) - Posb(1));
        Grad_AB_y = (Posa(3) - Posb(3)) / (Posa(2) - Posb(2));
        c_xz = Posa(3) - Grad_AB_z*Posa(1);
        c_yz = Posa(3) - Grad_AB_y*Posa(2);
        z_b = Grad_AB_z*x_b + c_xz;
        y_b = (z_b - c_yz)/Grad_AB_y;
        Posb = [x_b,y_b,z_b,1];
    end
end
if Posa(2) == Posb(2)
    Posb(2) = Posb(2) + 0.000000000000001;
end

Lambda_AB_xz = sqrt((Posb(1) - Posa(1))^2 + (Posb(3) - Posa(3))^2);

B_hat = Size*1.5;
% Theta
if Lambda_AB_xz == 0
    theta_xz = 0;
else
    theta_xz = asin(B_hat/Lambda_AB_xz);
end


T = [1 0 0 -Posa(1);
     0 1 0 -Posa(2);
     0 0 1 -Posa(3);
     0 0 0     1  ];
 
T_back = [1 0 0 Posa(1);
          0 1 0 Posa(2);
          0 0 1 Posa(3);
          0 0 0    1  ]; 
Posa;
Posb;
A = (T * Posa')';
B = (T * Posb')';

omega = -atan(sqrt(B(1)^2+B(2)^2)/B(3));
kapa = -atan(B(1)/B(2));       

Rz = [cos(kapa) sin(kapa) 0 0;
     -sin(kapa) cos(kapa) 0 0;
          0        0      1 0;
          0        0      0 1];
if Posa(2) <= Posb(2)
    Rx = [1      0           0      0;
          0  cos(omega)  sin(omega) 0;
          0 -sin(omega)  cos(omega) 0;
          0      0           0      1];
else
    Rx = [1      0            0     0;
          0  cos(omega) -sin(omega) 0;
          0  sin(omega)  cos(omega) 0;
          0      0           0      1];
end

R = Rx*Rz;
R_back = Rz'*Rx';
if B(:,1:2) ~= [0,0]
    B_new = (R* B');
else
    B_new = B;
end

Length = ceil(B_new(3));
if Length > 5 
    Length = 5;
elseif Length < -5
    Length = -5;
end

r = zeros(abs(Length/precision),1);
j = 1;
f = 1;
C = zeros(abs(Length/precision),3);
for i = 0:precision:abs(Length)
    r(j) = abs(Length - i) * sin(theta_xz); 
    
    % Centre of each disc
    Cx = 0;
    Cy = 0;
    if Length < 0 
        Cz = Length + i;
    else
        Cz = Length - i;
    end
    C(j,:) = [Cx,Cy,Cz];
    
    for m = -r(j):precision:r(j)
        for n = 0:pi/8:2*pi           
            x = m*sin(n);
            y = m*cos(n);
            z = 0;
            
            px = (C(j,1) + x);
            py = (C(j,2) + y);
            pz = (C(j,3) + z);            
            V_Coll = [px,py,pz,2] - A;
            Vel_Collide(f,:) = round(V_Coll/precision)*precision;
            %plot3(V_Coll(1)+A(1),V_Coll(2)+A(2),V_Coll(3)+A(3),'.k'); hold on
            f = f + 1;
        end        
    end                
    j = j + 1;
end



Vel_Collide_union = union(prev_collisionVelocity,Vel_Collide,'rows');
% Vel_Collide_union = unique(Vel_Collide_union,'rows');

if length(velocityA) < 4
    velocityA = [velocityA,1];
end

if velocityA(1:2) == [0 0]  
    vel_beg = velocityA;    
else
    vel_beg = velocityA;
    velocityA = (Rx * velocityA')';
    velocityA = (Rz * velocityA')';
    velocityA = round(velocityA/precision)*precision;
end

if currentVelocity(1:2) == [0 0]
    initialVelocity = currentVelocity;
else    
    initialVelocity = (Rz * [currentVelocity,1]');
    initialVelocity = (Rx * [currentVelocity,1]')';
    initialVelocity = round(initialVelocity/precision)*precision;
end

B_New = B_new';

% plot3(VelA(1)+A(1),VelA(2)+A(2),VelA(3)+A(3),'.b'); hold on
% plot3(PosA(1),PosA(2),PosA(3),'sr');
% plot3(PosB(1),PosB(2),PosB(3),'sb');
% plot3(Vel_start(1)+A(1),Vel_start(2)+A(2),Vel_start(3)+A(3),'.r');
% Cone(A(:,1:3),B_New(:,1:3),[0 B_hat],20,'m',0,0); hold on
% alpha(0.1)
velStart = 0;
max_vel = 5;
% If VelA a collision velocity calculate new velocity
if ismember(initialVelocity(1:3),Vel_Collide_union(:,1:3),'rows') == 0    
    velocityA = [currentVelocity,1];    
    velStart = 1;
elseif ismember(velocityA(1:3),Vel_Collide_union(:,1:3),'rows') == 1    
    A_tilda = A(:,1:3) + velocityA(:,1:3); 
    qx = A_tilda(1) - B_hat/2;
    qy = A_tilda(2) - B_hat/2;
    qz = A_tilda(3) - B_hat/2;
    Vel = zeros(29791,4);
    f = 1;
    
    for i = 0:precision:B_hat
       for j = 0:precision:B_hat
          for k = 0:precision:B_hat

             Q = [qx,qy,qz];
             if (Q - A(:,1:3)) <= (sqrt(max_vel^2 + max_vel^2))
                V = Q - A(:,1:3);
             end             
             Vel(f,1:3) = round(V/precision)*precision;
             if ismember(Vel(f,1:3),Vel_Collide_union(:,1:3),'rows') == 0
                 Vel_No_Collide(f,:) = round(Vel(f,:)/precision)*precision;
                 %plot3(A(1)+Vel(f,1),A(2)+Vel(f,2),A(3)+Vel(f,3),'.k'); hold on
             end
             f = f + 1;
             qx = qx + precision;
          end
          qx = A_tilda(1) - B_hat/2;
          qy = qy + precision;
       end
       qy = A_tilda(2) - B_hat/2;
       qz = qz + precision;
    end
    
    
    Near_dist = zeros(length(Vel_No_Collide),1);

    for m = 1:length(Vel_No_Collide)
        temp = velocityA - Vel_No_Collide(m,:);
        Temp = sqrt(temp(1)^2 + temp(2)^2 + temp(3)^2);
        Near_dist(m,:) = Temp;
    end
    [Closest,I] = min(Near_dist);
    New_Vel = Vel_No_Collide(I,:);
    velocityA = New_Vel;
end
if vel_beg(1:2) == [0 0]    
else
    velocityA = (Rz' * velocityA')';
    velocityA = (Rx' * velocityA')';
    velocityA = round(velocityA/precision)*precision;
end
if velStart == 1
    velocityA = (Rz' * velocityA')';
    velocityA = (Rx' * velocityA')';
    velocityA = round(velocityA/precision)*precision;
end

%%% Plot Collision Cone
% [X,Y,Z] = sphere;
% plot3(PosA(1),PosA(2),PosA(3),'sr'); hold on
% plot3(VelA(1)+PosA(1),VelA(2)+PosA(2),VelA(3)+PosA(3),'.b');
% plot3(Vel_Start(1)+PosA(1),Vel_Start(2)+PosA(2),Vel_Start(3)+PosA(3),'.r');
% Cone(Posa(:,1:3),Posb(:,1:3),[0 B_hat],20,colour,0,0); hold on
% alpha(0.1)
% surf(X+PosB(1),Y+PosB(2),Z+PosB(3));
% 
% xlim([-20 20]);
% ylim([-20 20]);
% zlim([-20 20]);
% grid on

% title('3D Velocity Obstacle')
% xlabel('x axis')
% ylabel('y axis')
% zlabel('z axis')

end
