function Agent_Avoid(Packet,Size)
%% Agent Avoid uses the Velocity Obstacle method
colour = [0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1;0   0   0.5;
0   0.5 0;0   0.5 0.5;0.5 0   0;0.5 0   0.5;0.5 0.5 0;0.5 0.5 0.5;
0   1   0.5;0   0.5 1;1   0   0.5;1   0.5 0;1   0.5 0.5;1   1   0.5;
1   0.5 1;0.5 0   1;0.5 1   0;0.5 1   0.5;0.5 0.5 1;0.5 1   1;0.7 0   0;
0   0.7 0;0   0   0.7;0   0.7 0.7;0.7 0   0.7;0.7 0.7 0.7];

    Vel_start = Packet(:,4:6);
    New_Vel = zeros(size(Packet,1),3);
    
    for runTime = 0:0
        
        for ID = 1:size(Packet,1)              
            New_Vel(ID,:) = Agent(ID,Packet,Vel_start,Size,colour(ID,:));                                                                       
        end
        Packet(:,1:3) = Packet(:,1:3) + New_Vel;
        Packet(:,4:6) = New_Vel;
    end      
end
