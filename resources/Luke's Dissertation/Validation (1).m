%% Validation Tests 2D
% Random Positions for Agent A and B
% Random Velocities for Agent A and B
clc 
cla
clear 

Size = 1;
LineColour = ['m','r','b','c','g','k'];

scenarios = 10000;
runTime = 20;
A = zeros([scenarios,6*runTime]);
B = zeros([scenarios,6*runTime]);
Collision = zeros([scenarios,runTime]);

for ID = 1:scenarios
    ID
    PosA = [randi([-10, 10],1,2),0];
    PosB = [randi([-10, 10],1,2),0];

    VelA = [randi([-5, 5],1,2),0];
    VelB = [randi([-5, 5],1,2),0];

    Vel_start = VelA;
    Vel_Collide_Prev = [0,0];
    
    %Log Position and Velocities
    Rnd = randi([1,6],1);
    
    for j = 0:runTime
        A(ID,j*6+1) = PosA(1);
        A(ID,j*6+2) = PosA(2);
        A(ID,j*6+3) = PosA(3);
        A(ID,j*6+4) = VelA(1);
        A(ID,j*6+5) = VelA(2);
        A(ID,j*6+6) = VelA(3);
        
        B(ID,j*6+1) = PosB(1);
        B(ID,j*6+2) = PosB(2);
        B(ID,j*6+3) = PosB(3);
        B(ID,j*6+4) = VelB(1);
        B(ID,j*6+5) = VelB(2);
        B(ID,j*6+6) = VelB(3);
        % Call Avoid Algorithm
        [VelA,Vel_Collide_Prev,Vel_No_Collide] = Avoid_Algorithm_2D(PosA,PosB,VelA,VelB,Size,Vel_Collide_Prev,Vel_start,LineColour(Rnd));
        PosA = PosA + VelA;
        PosB = PosB + VelB;
        if sqrt((PosA(1)-PosB(1))^2 + (PosA(2)-PosB(2))^2 + (PosA(3)-PosB(3))^2)  <= 1
           Collision(ID,j+1) = 1; 
        end
        [row, column] = find(Collision);
    end
end
%% Validation Tests 3D
% Random Positions for Agent A and B
% Random Velocities for Agent A and B
clc 
cla
clear 

Size = 1;
LineColour = ['m','r','b','c','g','k'];

scenarios = 50;
runTime = 20;
A = zeros([scenarios,6*runTime]);
B = zeros([scenarios,6*runTime]);
Collision = zeros([scenarios,runTime]);

for ID = 1:scenarios
    ID
    PosA = randi([-10, 10],1,3);
    PosB = randi([-10, 10],1,3);

    VelA = randi([-5, 5],1,3);
    VelB = randi([-5, 5],1,3);

    Vel_start = VelA;
    Vel_Collide_Prev = [0,0,0];
    
    %Log Position and Velocities
    Rnd = randi([1,6],1);
    
    for j = 0:runTime
        
        A(ID,j*6+1) = PosA(1);
        A(ID,j*6+2) = PosA(2);
        A(ID,j*6+3) = PosA(3);
        A(ID,j*6+4) = VelA(1);
        A(ID,j*6+5) = VelA(2);
        A(ID,j*6+6) = VelA(3);
        
        B(ID,j*6+1) = PosB(1);
        B(ID,j*6+2) = PosB(2);
        B(ID,j*6+3) = PosB(3);
        B(ID,j*6+4) = VelB(1);
        B(ID,j*6+5) = VelB(2);
        B(ID,j*6+6) = VelB(3);
        % Call Avoid Algorithm
        [VelA,Vel_Collide_Prev] = Avoid_Algorithm_3D_rotate(PosA,PosB,VelA,VelB,Size,Vel_Collide_Prev,Vel_start,LineColour(Rnd));      
        if sqrt((PosA(1)-PosB(1))^2 + (PosA(2)-PosB(2))^2 + (PosA(3)-PosB(3))^2)  <= 1
           Collision(ID,j+1) = 1; 
        end
        [row, column] = find(Collision);
        PosA = PosA + VelA;
        PosB = PosB + VelB;
    end
end

%% Collision Course
% Circle crossover
%profile on
clc 
cla
clear 

Size = 1;
colour = [0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;0 0 0;0   0   0.5;
0   0.5 0;0   0.5 0.5;0.5 0   0;0.5 0   0.5;0.5 0.5 0;0.5 0.5 0.5;
0   1   0.5;0   0.5 1;1   0   0.5;1   0.5 0;1   0.5 0.5;1   1   0.5;
1   0.5 1;0.5 0   1;0.5 1   0;0.5 1   0.5;0.5 0.5 1;0.5 1   1;0.7 0   0;
0   0.7 0;0   0   0.7;0   0.7 0.7;0.7 0   0.7;0.7 0.7 0.7];
Rnd = randi([1 23],1);

%%% Number of Agents %%%
n = 4;


phi = 2*pi / n;
r = 20;
m = 1;
Pos = zeros(n,3);
Vel = zeros(n,3);
Packet = zeros(n,6);

precision = 0.5;

for k = 0:phi:2*pi-phi
    x = r*sin(k);
    y = r*cos(k);
    u = -x/10;
    v = -y/10;
    Pos(m,:) = round([x,y,0]/precision)*precision;
    Vel(m,:) = round([u,v,0]/precision)*precision;
    Packet(m,:) = [Pos(m,:),Vel(m,:)];
    m = m + 1;
end

Vel_start = Vel; 
tests = 1;

runTime = 20*10;
Collision = [];
posOld = zeros(n,runTime*3);
for t = 1:tests
    t
    for r = 0:runTime-1
        r
        cla
        for ID = 1:n 
            
            [New_Vel(ID,:),Collision] = Agent(ID,Packet,Vel_start,Size,colour(ID,:),precision,Collision);
           
            %FinalVel = plot(New_Vel(ID,1)+Packet(ID,1),New_Vel(ID,2)+Packet(ID,2),'.g');
            posOld(ID,((r*3)+1):((r*3)+3)) = Packet(ID,1:3);
            
            plot(Packet(ID,1),Packet(ID,2),'Marker','o','MarkerFaceColor',colour(ID,:),'MarkerEdgeColor',colour(ID,:)); hold on
            for k = 0:r-1
                line([posOld(ID,(k*3)+1),posOld(ID,(k*3)+4)],[posOld(ID,(k*3)+2),posOld(ID,(k*3)+5)],'color',colour(ID,:));hold on                
            end
%             quiver(Packet(ID,1),Packet(ID,2),New_Vel(ID,1),New_Vel(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','k');hold on
%             quiver(Packet(ID,1),Packet(ID,2),Vel_start(ID,1),Vel_start(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','r');

        end 
        New_Vel;
        Packet(:,1:3) = Packet(:,1:3) + 0.1*New_Vel;
        Packet(:,4:6) = New_Vel;
                 
        grid on
        xlim([-25 25])
        ylim([-25 25])
        str = sprintf('%d Agents on a Collison Course',n);
        title(str)
        xlabel('x-axis')
        ylabel('y-axis')
        
%         filename = '4_Agents_0.1.gif';
%         drawnow
%         frame = getframe(1);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if r == 0
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.15);
%         end
%         pause
    end
end
%profile viewer
%profile off
%% Collision Course 3D
profile on
clc 
cla
clear 

Size = 1;
colour = [0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1;0   0   0.5;
0   0.5 0;0   0.5 0.5;0.5 0   0;0.5 0   0.5;0.5 0.5 0;0.5 0.5 0.5;
0   1   0.5;0   0.5 1;1   0   0.5;1   0.5 0;1   0.5 0.5;1   1   0.5;
1   0.5 1;0.5 0   1;0.5 1   0;0.5 1   0.5;0.5 0.5 1;0.5 1   1;0.7 0   0;
0   0.7 0;0   0   0.7;0   0.7 0.7;0.7 0   0.7;0.7 0.7 0.7];

%%% Number of Agents %%%
n = 4;
precision = 0.1;

phi = 2*pi/n;
r = 20;
m = 1;
Pos = zeros(n,3);
Vel = zeros(n,3);
Packet = zeros(n,6);

rnd = randi([1,6],1);

for k = 0:phi:2*pi-phi
    x = r*sin(k);
    y = r*cos(k);
    z = 5;
    u = -x/10;
    v = -y/10;
    w = 0;
    Pos(m,:) = round([x,y,z]/precision)*precision;
    Vel(m,:) = round([u,v,w]/precision)*precision;
    Packet(m,:) = [Pos(m,:),Vel(m,:)];    
    m = m + 1;
end
    
    
    

Vel_start = Vel; 
tests = 1;

runTime = 20;
Collision = [];
posOld = zeros(n,runTime*3);
for t = 1:tests
    t
    for r = 0:runTime-1
        r
        cla
        for ID = 1:n 
            
            New_Vel(ID,:) = Agent(ID,Packet,Vel_start,Size,colour(ID,:),precision,Collision);
            posOld(ID,((r*3)+1):((r*3)+3)) = Packet(ID,1:3);
            
            plot3(Packet(ID,1),Packet(ID,2),Packet(ID,3),'Marker','o','MarkerFaceColor',colour(ID,:),'MarkerEdgeColor',colour(ID,:)); hold on
            for k = 0:r-1
                plot3(posOld(ID,(k*3)+1),posOld(ID,(k*3)+2),posOld(ID,(k*3)+3),'LineStyle','-','Color',colour(ID,:))
            end
            
            %plot3(New_Vel(ID,1)+Packet(ID,1),New_Vel(ID,2)+Packet(ID,2),New_Vel(ID,3)+Packet(ID,3),'.g');
            
            
            
%             quiver(Packet(ID,1),Packet(ID,2),New_Vel(ID,1),New_Vel(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','k');hold on
%             quiver(Packet(ID,1),Packet(ID,2),Vel_start(ID,1),Vel_start(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','r');
            
        end                    
        Packet(:,1:3) = Packet(:,1:3) + 0.1*New_Vel(:,1:3);
        Packet(:,4:6) = New_Vel(:,1:3); 
         
        filename = '3D_crossover.gif';
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if r == 0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.15);
        end
        
        xlim([-25 25]);
        ylim([-25 25]);
        zlim([0 10]);
        grid on        
        title('3D Velocity Obstacle')
        xlabel('x axis')
        ylabel('y axis')
        zlabel('z axis')
        pause
    end
end
profile viewer
profile off
%% Collision Course
% Line crossover
%profile on
clc 
cla
clear 

Size = 1;
colour = [0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;0 0 0;0   0   0.5;
0   0.5 0;0   0.5 0.5;0.5 0   0;0.5 0   0.5;0.5 0.5 0;0.5 0.5 0.5;
0   1   0.5;0   0.5 1;1   0   0.5;1   0.5 0;1   0.5 0.5;1   1   0.5;
1   0.5 1;0.5 0   1;0.5 1   0;0.5 1   0.5;0.5 0.5 1;0.5 1   1;0.7 0   0;
0   0.7 0;0   0   0.7;0   0.7 0.7;0.7 0   0.7;0.7 0.7 0.7];
Rnd = randi([1 23],1);

%%% Number of Agents %%%
n = 8;


phi = 2*pi / n;
r = 20;
m = 1;
Pos = zeros(n,3);
Vel = zeros(n,3);
Packet = zeros(n,6);

precision = 0.5;

Pos = [
  0 -10 0;
-10  -6 0;
10  -4 0;
-10  -2 0;
10   0 0;
-10   2 0;
10   4 0;
-10   6 0];
Vel = [
  0  2  0;
  3  0  0;
 -2.8  0  0;
 2.6  0  0;
 -2.4  0  0;
 2.2  0  0;
 -2.0  0  0;
 1.8  0  0]; 

for a = 1:n
    Packet(a,:) = [Pos(a,:),Vel(a,:)];
end
Vel_start = Vel; 
tests = 1;

runTime = 20*10;
Collision = [];
posOld = zeros(n,runTime*3);
for t = 1:tests
    t
    for r = 0:runTime-1
        r
        cla
        for ID = 1:n 
            
            [New_Vel(ID,:),Collision] = Agent(ID,Packet,Vel_start,Size,colour(ID,:),precision,Collision);
           
            %FinalVel = plot(New_Vel(ID,1)+Packet(ID,1),New_Vel(ID,2)+Packet(ID,2),'.g');
            posOld(ID,((r*3)+1):((r*3)+3)) = Packet(ID,1:3);
            
            plot(Packet(ID,1),Packet(ID,2),'Marker','o','MarkerFaceColor',colour(ID,:),'MarkerEdgeColor',colour(ID,:)); hold on
            for k = 0:r-1
                line([posOld(ID,(k*3)+1),posOld(ID,(k*3)+4)],[posOld(ID,(k*3)+2),posOld(ID,(k*3)+5)],'color',colour(ID,:));hold on                
            end
%             quiver(Packet(ID,1),Packet(ID,2),New_Vel(ID,1),New_Vel(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','k');hold on
%             quiver(Packet(ID,1),Packet(ID,2),Vel_start(ID,1),Vel_start(ID,2),'MaxHeadSize',5,'LineWidth',2,'Color','r');

        end 
        New_Vel;
        Packet(:,1:3) = Packet(:,1:3) + 0.1*New_Vel;
        Packet(:,4:6) = New_Vel;
                 
        grid on
        xlim([-10 10])
        ylim([-10 10])
        str = sprintf('%d Agents Cross Over',n);
        title(str)
        xlabel('x-axis')
        ylabel('y-axis')
        
        filename = '8_Agents_crossover.gif';
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if r == 0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.15);
        end
%         pause
    end
end
%profile viewer
%profile off