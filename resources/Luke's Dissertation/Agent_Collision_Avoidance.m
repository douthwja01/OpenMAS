%% Collision Avoidance
profile on

clear 
cla
clc
addpath('agents');  
state  = [0;0;0;0;0;0];

state1 = [0;0;-2; 0;1;2];
state2 = [5;0;3; -2;0;0];
state3 = [5;-5;1.5; -2;2;0]; 


ID = zeros(1,3);
for i=1:length(ID)
    agent(i) = Basic_Drone(''); %ID{i} 
    agent(i) = agent(i).initialise(state);
end

agent(1) = agent(1).update_state(state1);
agent(2) = agent(2).update_state(state2);
agent(3) = agent(3).update_state(state3);
 
Packet = zeros(length(ID),6);

for i=1:length(ID)
    Pos = agent(i).position;
    Vel = agent(i).velocity;
    Size = agent(i).Size;
    Packet(i,:) = [Pos',Vel'];
end

Agent_Avoid(Packet,Size);
profile viewer
profile off