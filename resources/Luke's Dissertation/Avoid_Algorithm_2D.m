function [VelA,Vel_Collide_union] = Avoid_Algorithm_2D(PosA,PosB,VelA,VelB,Size,Vel_Collide_Prev,Vel_start,LineColour,precision)
%% Avoid Algorithm 2D
% Calculate Distance between 2 agents
Lambda_AB = sqrt((PosB(1) - PosA(1))^2+(PosB(2) - PosA(2))^2);  
% Calculate Radius of Agent 2, B_hat
B_hat = Size*1.5;
% Calculate length of tanget from point A to circle
Lambda_AC = sqrt(Lambda_AB^2 - B_hat^2);
Lambda_AD = Lambda_AC;
% Calculate Angles Phi, Theta, Gamma, Beta
% Phi
phi = acos(abs(PosA(1) - PosB(1))/Lambda_AB);
% Theta
theta = asin(B_hat/Lambda_AB);     
% Gamma
gamma = theta + phi;
% Beta
beta = gamma - 2*theta;
% Shift A by velocity of B to find new base point a
x_a = PosA(1) + VelB(1);
y_a = PosA(2) + VelB(2);
% Find Location of Agent B relative to Agent A
% Extend the tangents points past Agent B
% Quad 1
if (PosB(1) >= PosA(1)) && (PosB(2) >= PosA(2))
    if (PosB(1) < PosA(1)+B_hat)
        %ypos_xpos zone
        x_c = PosA(1) + (Lambda_AC * cos(gamma));
        y_c = PosA(2) + (Lambda_AC * sin(gamma));
        x_d = PosA(1) + (Lambda_AD * cos(beta));
        y_d = PosA(2) + (Lambda_AD * sin(beta)); 
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c - 3 + VelB(1);
        x_d = x_d + 3 + VelB(1);            
    else
        x_c = PosA(1) + (Lambda_AC * cos(gamma));
        y_c = PosA(2) + (Lambda_AC * sin(gamma));
        x_d = PosA(1) + (Lambda_AD * cos(beta));
        y_d = PosA(2) + (Lambda_AD * sin(beta));
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c + 3 + VelB(1);
        x_d = x_d + 3 + VelB(1);            
     end
% Quad 2    
elseif (PosB(1) >= PosA(1)) && (PosB(2) <= PosA(2))
    if (PosB(1) < PosA(1)+B_hat)
        %yneg_xpos
        x_c = PosA(1) + (Lambda_AC * cos(beta));
        y_c = PosA(2) - (Lambda_AC * sin(beta));
        x_d = PosA(1) + (Lambda_AD * cos(gamma));
        y_d = PosA(2) - (Lambda_AD * sin(gamma)); 
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c + 3 + VelB(1);
        x_d = x_d - 3 + VelB(1);            
    else
        x_c = PosA(1) + (Lambda_AC * cos(beta));
        y_c = PosA(2) - (Lambda_AC * sin(beta));
        x_d = PosA(1) + (Lambda_AD * cos(gamma));
        y_d = PosA(2) - (Lambda_AD * sin(gamma)); 
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c + 3 + VelB(1);
        x_d = x_d + 3 + VelB(1);            
    end
% Quad 3    
elseif (PosB(1) <= PosA(1)) && (PosB(2) <= PosA(2))
    if (PosB(1) > PosA(1)-B_hat)
        %yneg_xneg
        x_c = PosA(1) - (Lambda_AC * cos(gamma));
        y_c = PosA(2) - (Lambda_AC * sin(gamma));
        x_d = PosA(1) - (Lambda_AD * cos(beta));
        y_d = PosA(2) - (Lambda_AD * sin(beta));
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c + 3 + VelB(1);
        x_d = x_d - 3 + VelB(1);               
    else    
        x_c = PosA(1) - (Lambda_AC * cos(gamma));
        y_c = PosA(2) - (Lambda_AC * sin(gamma));
        x_d = PosA(1) - (Lambda_AD * cos(beta));
        y_d = PosA(2) - (Lambda_AD * sin(beta));
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c - 3 + VelB(1);
        x_d = x_d - 3 + VelB(1);            
    end
% Quad 4    
elseif (PosB(1) < PosA(1)) && (PosB(2) >= PosA(2))
    if (PosB(1) > PosA(1)-B_hat)
        % ypos_xneg
        x_c = PosA(1) - (Lambda_AC * cos(beta));
        y_c = PosA(2) + (Lambda_AC * sin(beta));
        x_d = PosA(1) - (Lambda_AD * cos(gamma));
        y_d = PosA(2) + (Lambda_AD * sin(gamma));
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c - 3 + VelB(1);
        x_d = x_d + 3 + VelB(1);              
    else    
        x_c = PosA(1) - (Lambda_AC * cos(beta));
        y_c = PosA(2) + (Lambda_AC * sin(beta));
        x_d = PosA(1) - (Lambda_AD * cos(gamma));
        y_d = PosA(2) + (Lambda_AD * sin(gamma));
        Grad_AC_xy = ((y_c+VelB(2)) - y_a)/((x_c+VelB(1)) - x_a);
        Grad_AD_xy = ((y_d+VelB(2)) - y_a)/((x_d+VelB(1)) - x_a);
        x_c = x_c - 3 + VelB(1);
        x_d = x_d - 3 + VelB(1);            
    end
end
% Extending the y coordinates of C and D by using equation of a line y = mx+x    
c_AC_xy = y_a - Grad_AC_xy*x_a;            
c_AD_xy = y_a - Grad_AD_xy*x_a;           
y_c = Grad_AC_xy*(x_c) + c_AC_xy;
y_d = Grad_AD_xy*(x_d) + c_AD_xy;
% Positions of Points Shifted by VelB
Posa = [PosA(1)+VelB(1),PosA(2)+VelB(2),PosA(3)+VelB(3)];
PosB = [PosB(1),PosB(2),PosB(3)];
PosC = [x_c,y_c];
PosD = [x_d,y_d];
% Set up scan starting at -5m away
x = PosA(1) - 5; 
y = PosA(2) - 5;                          
% Create size of variables                                          
Vel_Collide = ones(100,2)*10;
Vel_No_Collide = ones(500,2)*10;
a = 1;
b = 1;
p1 = Posa; 
p2 = PosC; 
p3 = PosD;
for i = 0:10/precision
    for j = 0:10/precision 
        p = [x,y];
        % Calculating Barycentric Coordinates
        alpha = ((p2(2) - p3(2))*(p(1) - p3(1)) + (p3(1) - p2(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)));
        eta  = ((p3(2) - p1(2))*(p(1) - p3(1)) + (p1(1) - p3(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)));
        zeta = 1 - alpha - eta;       
        if ((alpha >= 0) && (eta >= 0) && (zeta >= 0) && (alpha <= 1) && (eta <= 1) && (zeta <= 1))
            P_Col = [x,y];
            V_Col = P_Col - [PosA(1),PosA(2)];
            Vel_Collide(a,:) = round(V_Col/precision)*precision;        
            a = a + 1;            
        else 
            P_NoCol = [x,y];
            V_NoCol = P_NoCol - [PosA(1),PosA(2)];
            Vel_No_Collide(b,:) = V_NoCol;
            b = b + 1;            
        end 
        y = y + precision;
    end
    y = PosA(2) - 5;
    x = x + precision;
end
% Create union of past velocity obstacles
Vel_Collide_union = union(Vel_Collide_Prev,Vel_Collide,'rows');
% Set size of variables
Near_dist = ones(length(Vel_No_Collide),1)*10;
temp = zeros(length(Vel_No_Collide),2);
n = 1;
% Check if VelA is a Collision Velocity
while ismember(VelA(1:2),Vel_Collide_union,'rows') == 1
    if n > 1
        Vel_No_Collide(I,:) = [];
    end
    % Change Velocity to nearest value that isn't in the collision set   
    for m = 1:length(Vel_No_Collide)
        temp(m,:) = VelA(:,1:2) - Vel_No_Collide(m,:);
        Temp = sqrt(temp(m,1)^2 + temp(m,2)^2);
        Near_dist(m) = Temp;        
    end
    [Closest,I] = min(Near_dist);
    New_Vel = Vel_No_Collide(I,:);
    VelA(:,1:2) = New_Vel;
    n = n + 1;
end
% If Vel start not a collision velocity use vel start as new velocity
if ismember(Vel_start(:,1:2),Vel_Collide_union,'rows') == 0
    VelA(:,1:2) = Vel_start(:,1:2);    
end
%%% Plotting functions
plot(PosA(1),PosA(2),'or'); hold on
viscircles(PosA(:,1:2),Size/2,'EdgeColor',LineColour,'LineWidth',0.8);
viscircles(PosB(:,1:2),Size/2,'EdgeColor',LineColour,'LineWidth',0.8);
NewVel = plot(VelA(1)+PosA(1),VelA(2)+PosA(2),'.b');
StartVel = plot(Vel_start(1)+PosA(1),Vel_start(2)+PosA(2),'.r');
plot([Posa(1),PosC(1)],[Posa(2),PosC(2)],'Color',LineColour);
plot([Posa(1),PosD(1)],[Posa(2),PosD(2)],'Color',LineColour);
drawnow
refresh
xlim([-20 20]);
ylim([-20 20]);
grid on
title('Testing New Velocity Vector')
xlabel('x axis')
ylabel('y axis')
legend([StartVel,NewVel],' Start Velocity ','  New Velocity')
end