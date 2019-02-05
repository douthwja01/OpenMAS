%% Receives Packet of information about all Agents
%  Removes own information
%  Uses rest of packet on avoid algorithm
function [Vel,Collision] = Agent(ID,Packet,Vel_start,Size,colour,precision,Collision)
    
    % Packet(i,:) = [Pos',Vel'];

    Pos = Packet(ID,1:3);
    Vel = round(Packet(ID,4:6)/precision)*precision;
    Vel_start = round(Vel_start(ID,:)/precision)*precision;
    Packet(ID,:) = []; 
    
    
    Vel_Collide_Prev = [10,10];    
    Vel_Collide_Prev_3D = [10,10,10,10];
    
    
%     ID
    for n = 1:size(Packet,1)
%         n  
        [Vel,Vel_Collide_Prev] = Avoid_Algorithm_2D(Pos,Packet(n,1:3),Vel,Packet(n,4:6),Size,Vel_Collide_Prev,Vel_start,colour,precision);
%         [Vel,Vel_Collide_Prev_3D] = Avoid_Algorithm_3D_allRotate(Pos,Packet(n,1:3),Vel,Packet(n,4:6),Size,Vel_Collide_Prev_3D,Vel_start,colour,precision);
        
        if sqrt((Pos(1)-Packet(n,1))^2 + (Pos(2)-Packet(n,2))^2 + (Pos(3)-Packet(n,3))^2) <= Size
            Collision = [Collision;1];
        end
        
    end
    
end


% [Vel_1,Vel_Collide_Prev1_2D] = Avoid_Algorithm_2D(Pos_1,Pos_2,Vel_1,Vel_2,Size,Vel_Collide_Prev1_2D,Vel_1_start,'m');
% [Vel_1,Vel_Collide_Prev1_3D] = Avoid_Algorithm_3D(Pos_1,Pos_3,Vel_1,Vel_3,Size,Vel_Collide_Prev1_3D,Vel_1_start,'m');
% [Vel_1,Vel_Collide_Prev1_3D] = Avoid_Algorithm_3D_rotate(Pos_1,Pos_2,Vel_1,Vel_2,Size,Vel_Collide_Prev1_3D,Vel_1_start,'m');        
% [Vel_1,Vel_Collide_Prev1_3D] = Avoid_Algorithm_3D_allRotate(Pos_1,Pos_2,Vel_1,Vel_2,Size,Vel_Collide_Prev1_3D,Vel_1_start,'m');        