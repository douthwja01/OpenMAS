%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quiver test
PosA = [0,0,0];
PosB = [0,10,0];
VelB = [0,0,0];
Posa = PosA;
cla
% Calculate Distance between 2 agents
Lambda_AB_xy = sqrt((PosB(1) - PosA(1))^2+(PosB(2) - PosA(2))^2);
Lambda_AB_xz = sqrt((PosB(1) - PosA(1))^2+(PosB(3) - PosA(3))^2);
% Calculate Radius of Agent 2, B_hat
B_hat = 1.5;
% Calculate length of tanget from point A to circle
Lambda_AC_xy = sqrt(Lambda_AB_xy^2 - B_hat^2);
Lambda_AC_xz = sqrt(Lambda_AB_xz^2 - B_hat^2);
Lambda_AD_xy = Lambda_AC_xy;
Lambda_AD_xz = Lambda_AC_xz;
Lambda_AE_xz = sqrt(Lambda_AB_xz^2 - B_hat^2);
Lambda_AF_xz = Lambda_AE_xz;

% Calculate Angles Phi, Theta, Gamma, Beta
% Phi
phi_xy = acos(abs(PosA(1) - PosB(1))/Lambda_AB_xy);
phi_xz = acos(abs(PosA(1) - PosB(1))/Lambda_AB_xz);
% Theta
theta_xy = asin(B_hat/Lambda_AB_xy);     
theta_xz = asin(B_hat/Lambda_AB_xz);
% Gamma
gamma_xy = theta_xy + phi_xy;
gamma_xz = theta_xz + phi_xz;
% Beta
beta_xy = gamma_xy - 2*theta_xy;
beta_xz = gamma_xz - 2*theta_xz;

x_a = PosA(1);
y_a = PosA(2);
z_a = PosA(3);
Quad1 = 0;
Quad2 = 0;
Quad3 = 0;
Quad4 = 0;
ypos = 0;
yneg = 0;
xpos = 0;
xneg = 0;
%Location of Agent B
% Quad 1
if (PosB(1) > PosA(1)) && (PosB(2) > PosA(2))
    x_c = PosA(1) + (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) + (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) + (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) + (Lambda_AD_xy * sin(beta_xy)); 
    z_d = PosB(3);   
    Quad1 = 1;
% Quad 2    
elseif (PosB(1) > PosA(1)) && (PosB(2) < PosA(2))
    x_c = PosA(1) + (Lambda_AC_xy * cos(beta_xy));
    y_c = PosA(2) - (Lambda_AC_xy * sin(beta_xy));
    z_c = PosB(3);
    x_d = PosA(1) + (Lambda_AD_xy * cos(gamma_xy));
    y_d = PosA(2) - (Lambda_AD_xy * sin(gamma_xy)); 
    z_d = PosB(3);
    Quad2 = 1;
% Quad 3    
elseif (PosB(1) < PosA(1)) && (PosB(2) < PosA(2))
    x_c = PosA(1) - (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) - (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) - (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) - (Lambda_AD_xy * sin(beta_xy));
    z_d = PosB(3);
    Quad3 = 1;
% Quad 4    
elseif (PosB(1) < PosA(1)) && (PosB(2) > PosA(2))
    x_c = PosA(:,1) - (Lambda_AC_xy * cos(beta_xy));
    y_c = PosA(:,2) + (Lambda_AC_xy * sin(beta_xy));
    z_c = PosB(:,3);
    x_d = PosA(:,1) - (Lambda_AD_xy * cos(gamma_xy));
    y_d = PosA(:,2) + (Lambda_AD_xy * sin(gamma_xy));
    z_d = PosB(:,3);
    Quad4 = 1;
% y axis pos
elseif (PosB(1) == PosA(1)) && (PosB(2) > PosA(2))     
    x_c = PosA(1) - (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) + (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) + (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) + (Lambda_AD_xy * sin(beta_xy)); 
    z_d = PosB(3);
    ypos = 1;
% y axis neg
elseif (PosB(1) == PosA(1)) && (PosB(2) < PosA(2))     
    x_c = PosA(1) + (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) - (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) + (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) - (Lambda_AD_xy * sin(beta_xy)); 
    z_d = PosB(3);
    yneg = 1;
% x axis pos
elseif (PosB(1) > PosA(1)) && (PosB(2) == PosA(2))     
    x_c = PosA(1) + (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) + (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) + (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) - (Lambda_AD_xy * sin(beta_xy)); 
    z_d = PosB(3); 
    xpos = 1;
% x axis neg
elseif (PosB(1) < PosA(1)) && (PosB(2) == PosA(2))     
    x_c = PosA(1) - (Lambda_AC_xy * cos(gamma_xy));
    y_c = PosA(2) - (Lambda_AC_xy * sin(gamma_xy));
    z_c = PosB(3);
    x_d = PosA(1) - (Lambda_AD_xy * cos(beta_xy));
    y_d = PosA(2) + (Lambda_AD_xy * sin(beta_xy)); 
    z_d = PosB(3); 
    xneg = 1;
end

Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
Grad_AC_xz = ((z_c+VelB(3)) - z_a)/((x_c+VelB(1)) - x_a);
Grad_AD_xz = ((z_d+VelB(3)) - z_a)/((x_d+VelB(1)) - x_a);
c_AC_xy = y_a - Grad_AC_xy*x_a;            
c_AD_xy = y_a - Grad_AD_xy*x_a;
c_AC_xz = z_a - Grad_AC_xz*x_a;
c_AD_xz = z_a - Grad_AD_xz*x_a;
x_c = x_c + 3 + VelB(1,1);
x_d = x_d + 3 + VelB(1,1);            
y_c = Grad_AC_xy*(x_c) + c_AC_xy;
y_d = Grad_AD_xy*(x_d) + c_AD_xy;  
z_c = Grad_AC_xz*(x_c) + c_AC_xz;
z_d = Grad_AD_xz*(x_d) + c_AD_xz;

PosC = [x_c,y_c,z_c];
PosD = [x_d,y_d,z_d];
Lambda_AC_xy = sqrt((PosC(1) + Posa(1))^2 + (PosC(2) + Posa(2))^2);
Lambda_AC_xz = sqrt((PosC(1) + Posa(1))^2 + (PosC(3) + Posa(3))^2);
Lambda_AD_xy = sqrt((PosD(1) + Posa(1))^2 + (PosD(2) + Posa(2))^2);
Lambda_AD_xz = sqrt((PosD(1) + Posa(1))^2 + (PosD(3) + Posa(3))^2);

if Quad1 == 1
   alpha_xy = acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 1
elseif Quad2 == 1
   alpha_xy = pi - acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = pi - acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = pi - acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = pi - acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 2
elseif Quad3 == 1
   alpha_xy = pi + acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = pi + acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = pi + acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = pi + acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 3
elseif Quad4 == 1
   alpha_xy = 2*pi - acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = 2*pi - acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = 2*pi - acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = 2*pi - acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 4
elseif ypos == 1
   alpha_xy = -acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = -acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 'ypos'
elseif yneg == 1
   alpha_xy = acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 'yneg'  
elseif xpos == 1
   alpha_xy = acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = pi - acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = pi - acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 'xpos'   
elseif xneg == 1
   alpha_xy = pi + acos(abs(PosC(2) - Posa(2))/Lambda_AC_xy);
   alpha_xz = pi + acos(abs(PosC(3) - Posa(3))/Lambda_AC_xz);
   delta_xy = 2*pi - acos(abs(PosD(2) - Posa(2))/Lambda_AD_xy);
   delta_xz = 2*pi - acos(abs(PosD(3) - Posa(3))/Lambda_AD_xz);
   quad = 'xneg'   
end

x = PosA(1) - 5; 
y = PosA(2) - 5;                          % Sets area of 10.5m^2
z = PosA(3) - 5;                          % for velocity vectors

Vel_Collide = [0,0,0];
Vel_No_Collide = [0,0,0];
V_All = [0,0,0];


%for i = 0:20
%    for j = 0:20 
%        for k = 0:20
            %p = [x,y,z];
            p = [4 0 0;0 -4 0;-4 0 0;0 4 0];
            for i = 1:4
            lambda_AP_xy = sqrt((Posa(1) - p(i,1))^2 + (Posa(2) - p(i,2))^2);
            lambda_AP_xz = sqrt((Posa(1) - p(i,1))^2 + (Posa(3) - p(i,3))^2);
            
            if (p(i,1) >= PosA(1)) && (p(i,2) >= PosA(2))
               zeta_xy = acos(abs(p(i,2) - Posa(2)) / lambda_AP_xy);
               zeta_xz = acos(abs(p(i,3) - Posa(3)) / lambda_AP_xz);
               p_quad = 1
            elseif (p(i,1) >= PosA(1)) && (p(i,2) < PosA(2))
               zeta_xy = pi - acos(abs(p(i,2) - Posa(2)) / lambda_AP_xy);
               zeta_xz = pi - acos(abs(p(i,3) - Posa(3)) / lambda_AP_xz);
               p_quad = 2
            elseif (p(i,1) <= PosA(1)) && (p(i,2) <= PosA(2))
               zeta_xy = pi + acos(abs(p(i,2) - Posa(2)) / lambda_AP_xy);
               zeta_xz = pi + acos(abs(p(i,3) - Posa(3)) / lambda_AP_xz);
               p_quad = 3
            elseif (p(i,1) <= PosA(1)) && (p(i,2) > PosA(2))
               zeta_xy = 2*pi - acos(abs(p(i,2) - Posa(2)) / lambda_AP_xy);
               zeta_xz = 2*pi - acos(abs(p(i,3) - Posa(3)) / lambda_AP_xz);
               p_quad = 4
            end
            
            zeta_xy
            alpha_xy
            delta_xy
            zeta_xz
            alpha_xz
            delta_xz

            if (((zeta_xy >= alpha_xy) && (zeta_xy <= delta_xy)) && ((zeta_xz >= alpha_xz) && (zeta_xz <= delta_xz) ))
                P_Col = p(i,:)
                V_Col = P_Col - PosA;
                Vel_Collide = [Vel_Collide;V_Col];
                V_All = [V_All;V_Col];
                quiver3(PosA(1),PosA(2),PosA(3),V_Col(1),V_Col(2),V_Col(3),'Linewidth',1);
            else 
                P_NoCol = p(i,:);
                V_NoCol = P_NoCol - PosA;
                Vel_No_Collide = [Vel_No_Collide;V_NoCol];
                V_All = [V_All;V_NoCol];
                %quiver3(PosA(1),PosA(2),PosA(3),V_NoCol(1),V_NoCol(2),V_NoCol(3),'Linewidth',1);
            end          
            z = z + 0.5;
            end
         
%        end
        z = PosA(3) - 5;
        y = y + 0.5;
%    end
    y = PosA(2) - 5;
    x = x + 0.5;
%end
        
plot3(PosA(1),PosA(2),PosA(3),'rs'); hold on
Cone(PosA,PosB,[0 B_hat],20,'m',0,1);
Alpha(0.2)
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
grid on

%% Cone and sphere plot
A = [0 0 0];
B = [5 0 5];
r = [0 1];
[x,y,z] = sphere;

figure
Cone(A,B,r,20,'r',0,1); hold on
surf(x+5,y,z+5);
Alpha(0.9);
xlim([0 10])
ylim([-5 5])
zlim([0 10])
%% point in cone
Posa = [0,0,0];
PosA = [0,0,0];
PosB = [10,-2,0];
h = sqrt((PosA(1)-PosB(1))^2 + (PosA(2)-PosB(2))^2 + (PosA(3)-PosB(3))^2 );
r = 1;


x = PosA(1) - 5; 
y = PosA(2) - 5;                          % Sets area of 10.5m^2
z = PosA(3) - 5;                          % for velocity vectors

Vel_Collide = [];
Vel_No_Collide = [];
V_All = [];



for i = 0:10
    for j = 0:10 
        for k = 0:10
            p = [x,y,z];

            
            if (((zeta_xy >= beta_xy) && (zeta_xy <= gamma_xy)) && ((zeta_xz >= beta_xz) && (zeta_xz <= gamma_xz)))
                P_Col = [x,y,z];
                V_Col = P_Col - PosA;
                Vel_Collide = [Vel_Collide;V_Col];
                quiver3(PosA(1),PosA(2),PosA(3),V_Col(1),V_Col(2),V_Col(3),'Linewidth',1);
            else 
                P_NoCol = [x,y,z];
                V_NoCol = P_NoCol - PosA;
                Vel_No_Collide = [Vel_No_Collide;V_NoCol];
            end          
            z = z + 1;
        end
        z = PosA(3) - 5;
        y = y + 1;
    end
    y = PosA(2) - 5;
    x = x + 1;
end


plot3(PosA(1),PosA(2),PosA(3),'rs'); hold on
Cone(PosA,PosB,[0 r],20,'m',0,1);
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
grid on

%% Cone Function
clear
clc
cla
PosA = [0,0,0];
PosB = [-3,-3,10];
B_hat = 1;

Lambda_AB = sqrt((PosB(1) - PosA(1))^2+(PosB(2) - PosA(2))^2 + (PosB(3) - PosA(3))^2);
theta = asin(B_hat/Lambda_AB); 
Quad0_xz = 0;
Quad1_xz = 0;
Quad2_xz = 0;
Quad3_xz = 0;
Quad4_xz = 0;
Quad0_yz = 0;
Quad1_yz = 0;
Quad2_yz = 0;
Quad3_yz = 0;
Quad4_yz = 0;
%Location of Agent B xz plane
if (PosB(1) >= PosA(1)) && (PosB(3) >= PosA(3))   
    Quad1_xz = 1 
elseif (PosB(1) > PosA(1)) && (PosB(3) < PosA(3))
    Quad2_xz = 1
elseif (PosB(1) <= PosA(1)) && (PosB(3) < PosA(3))
    Quad3_xz = 1
elseif (PosB(1) < PosA(1)) && (PosB(3) >= PosA(3))
    Quad4_xz = 1
elseif (PosB(1) == PosA(1)) && (PosB(3) == PosA(3))
    Quad0_xz = 1
end
%Location of Agent B yz plane
if (PosB(2) >= PosA(2)) && (PosB(3) >= PosA(3))   
    Quad1_yz = 1 
elseif (PosB(2) > PosA(2)) && (PosB(3) < PosA(3))
    Quad2_yz = 1
elseif (PosB(2) <= PosA(2)) && (PosB(3) < PosA(3))
    Quad3_yz = 1
elseif (PosB(2) < PosA(2)) && (PosB(3) >= PosA(3))
    Quad4_yz = 1
elseif (PosB(2) == PosA(2)) && (PosB(3) == PosA(3))
    Quad0_yz = 1
end
Length = round(Lambda_AB);
r = zeros(Length,1);
j = 1;
p = [];
n = [];
for i = 0:1:Lambda_AB
    
    r(j) = (Lambda_AB - i) * sin(theta);   
     
    if Quad1_xz == 1
        psi_xz = asin(abs(PosB(1) - PosA(1))/Lambda_AB) ;
    elseif Quad2_xz == 1
        psi_xz = pi - asin(abs(PosB(1) - PosA(1))/Lambda_AB) ;
    elseif Quad3_xz == 1
        psi_xz = pi + asin(abs(PosB(1) - PosA(1))/Lambda_AB) ;
    elseif Quad4_xz == 1
        psi_xz = 2*pi - asin(abs(PosB(1) - PosA(1))/Lambda_AB) ;
    elseif Quad0_xz == 1
        psi_xz = 0;
    end
    
    if Quad1_yz == 1
        psi_yz = asin(abs(PosB(2) - PosA(2))/Lambda_AB) ;
    elseif Quad2_yz == 1
        psi_yz = pi - asin(abs(PosB(2) - PosA(2))/Lambda_AB) ;
    elseif Quad3_yz == 1
        psi_yz = pi + asin(abs(PosB(2) - PosA(2))/Lambda_AB) ;
    elseif Quad4_yz == 1
        psi_yz = 2*pi - asin(abs(PosB(2) - PosA(2))/Lambda_AB) ; 
    elseif Quad0_yz == 1
        psi_yz = 0;
    end
    
    % Centre of each disk
    
    Cx = (Lambda_AB - i)*sin(psi_xz);
    Cy = (Lambda_AB - i)*sin(psi_yz);
    Cz = (Lambda_AB - i)*cos(psi_xz);
    C = [Cx,Cy,Cz];
    radius = r(j);
    
    for m = -r(j):0.2:r(j)
        
        for n = 0:pi/4:2*pi   
            x = m*sin(n);
            y = m*cos(n);
            z = 0;
            
            px = Cx + x;
            py = Cy + y;
            pz = Cz + z;
            plot3(px,py,pz,'.k'); hold on
            p = [p;px,py,pz];
        end        
    end        
    plot3(Cx,Cy,Cz,'ok'); hold on
    j = j + 1;
end

Cone(PosA,PosB,[0 B_hat],20,'r',0,0);
Alpha(0.1);
xlim([-5 5])
ylim([-5 5])
zlim([0 12])

xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
title('Filling Cone with Collision Points ')
grid on


%%
cla
Posa = [1 0 0];
Posb = [0 0 0];
B_hat = 1.5;

plot3(Posa(1),Posa(2),Posa(3),'ob'); hold on
plot3(Posb(1),Posb(2),Posb(3),'sg')
Cone(Posa,Posb,[0 B_hat],20,'r',0,1);
Alpha(0.2)
grid on

%%  Velocity Obstacle Tests
clear 
cla
clc
Size = 1;

Pos_1 = [0,0,0];
Pos_2 = [5,5,0;
         5,0,0;
         0,5,0;
         5,-5,0;
        -5,5,0;
        -5,0,0;
         0,-5,0;
        -5,-5,0];
Vel_1 = [0,0,0];
Vel_2 = [1,0,0;
         1,1,0;
         0,1,0;
         0,-1,0;
        -1,-1,0;
        -1,0,0;
        -1,1,0;
         1,-1,0];

Vel_Collide_Prev1_2D    = [0,0];       
Vel_1_start = Vel_1;        
Colours = ['y','m','c','r','g','b','k','k'];
for i = 1:length(Pos_2)
%     [Vel_1,Vel_Collide_Prev1_2D] = Avoid_Algorithm_2D(Pos_1,Pos_2(i,:),Vel_1,Vel_2,Size,Vel_Collide_Prev1_2D,Vel_1_start,'m');
    [Vel_1,Vel_Collide_Prev1_2D] = Avoid_Algorithm_2D(Pos_1,Pos_2(2,:),Vel_1,Vel_2(i,:),Size,Vel_Collide_Prev1_2D,Vel_1_start,Colours(i));
end
Pos_1 = Pos_1 + Vel_1;
Pos_2 = Pos_2 + Vel_2;

% title('Testing Velocity Obstacle Orientation')
title('Testing Velocity Obstacle Start Point')
%% Transformation Matrices
cla
clf
clc
clear
Vela = [0,2,5,1];
Posa = [0,0,0,1];
Posb = [-5,0,0,1];
if Posa(2) == Posb(2)
    Posb(2) = Posb(2) + 0.000000000000001;
end
lambda_AB = sqrt((Posb(1) - Posa(1))^2 + (Posb(2) - Posa(2))^2 + (Posb(3) - Posa(3))^2);

Lambda_AB_xz = sqrt((Posb(1) - Posa(1))^2 + (Posb(3) - Posa(3))^2);
Vel_Collide = [0,0,0];
V_rot = [];
V_tran = [];
B_hat = 1.5;
% Theta
if Lambda_AB_xz == 0
    theta_xy = 0;
else
    theta_xy = asin(B_hat/Lambda_AB_xz);
end

T = [1 0 0 -Posa(1);
     0 1 0 -Posa(2);
     0 0 1 -Posa(3);
     0 0 0     1  ];
T_back = [1 0 0 Posa(1);
          0 1 0 Posa(2);
          0 0 1 Posa(3);
          0 0 0    1  ];    
     
A = (T * Posa')'
B = (T * Posb')'


Alpha = 0;
beta = acos(B(3)/lambda_AB);
gamma = atan(B(2)/B(1));

omega = -atan(sqrt(B(1)^2+B(2)^2)/B(3));
phi = 0;
kapa = -atan(B(1)/B(2));       

Mz = [cos(kapa) sin(kapa) 0 0;
     -sin(kapa) cos(kapa) 0 0;
          0        0      1 0;
          0        0      0 1];
if Posa(2) <= Posb(2)
    Mx = [1      0           0      0;
          0  cos(omega)  sin(omega) 0;
          0 -sin(omega)  cos(omega) 0;
          0      0           0      1];
else
    Mx = [1      0            0     0;
          0  cos(omega) -sin(omega) 0;
          0  sin(omega)  cos(omega) 0;
          0      0           0      1];
end 

M = Mx*Mz; 
M_back = Mz'*Mx';

B_new = (M* B');

Length = ceil(B_new(3));

r = zeros(Length,1);
j = 1;
acc = 0.5;
C = zeros(Length,3);
for i = 0:acc:Length    
    r(j) = (Length - i) * sin(theta_xy); 
    
    % Centre of each disc
    Cx = 0;
    Cy = 0;
    Cz = Length - i;
    C(j,:) = [Cx,Cy,Cz];
        
    for m = -r(j):acc:r(j)
        for n = 0:pi/8:2*pi   
            x = m*sin(n);
            y = m*cos(n);
            z = 0;
            
            px = (C(j,1) + x);
            py = (C(j,2) + y);
            pz = (C(j,3) + z);            
            V_Coll = [px,py,pz,2] - A;
        
            
%             plot3(V_Coll_trans(1),V_Coll_trans(2),V_Coll_trans(3),'.k'); hold on
%             plot3(V_Coll(1),V_Coll(2),V_Coll(3),'.m'); hold on
        end        
    end        
    
    
    j = j + 1;
end

% Rotate back and translate
B_rotback = M_back* B_new;
B_trans = T_back * B_rotback;





plot3(Posa(1),Posa(2),Posa(3),'or'); hold on
plot3(Posb(1),Posb(2),Posb(3),'ob')
plot3(A(1),A(2),A(3),'sr')
plot3(B(1),B(2),B(3),'sb')
plot3(B_new(1),B_new(2),B_new(3),'db')
plot3(A(1),A(2),A(3),'dr')
orig = [Posa;Posb];
shift = [A;B];
new = [A;B_new'];
orig = line(orig(:,1),orig(:,2),orig(:,3));
new = line(new(:,1),new(:,2),new(:,3),'Color','r');
shift = line(shift(:,1),shift(:,2),shift(:,3),'Color','g');
B_new = B_new';
% Cone(Posa(:,1:3),Posb(:,1:3),[0 B_hat],20,'k',0,0); hold on
% Cone(A(:,1:3),B_new(:,1:3),[0 B_hat],20,'m',0,0);
alpha(0.1)
grid on;
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
xlim([-5 5])
ylim([-5 5])
zlim([0 10])
%legend([orig,shift,new],'Original','Shifted','Rotated')
title('Collision Points Inside Cone')




