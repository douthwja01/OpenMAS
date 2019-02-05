% UPDATE AGENT STATE VECTOR [x;y;z;u;v;w;psi;the;phi;p;q;r]
function newState = stateDynamics_velocities(obj,dt,linearVelocity,angularVelocity)

% LOCAL STATE:
% [x; y; z; phi; theta; psi; u; v; w; phi_dot; theta_dot; psi_dot]

dX = [linearVelocity;angularVelocity;zeros(6,1)];

% INTEGRATION MATRIX
dtMatrix = dt*ones(size(obj.localState));

% DEFINE THE NEW STATE VECTOR
newState = obj.localState + dtMatrix*dX;                   % Access the current state and append difference
% AMMEND THE NEW VELOCITIES
newState(7:12) = [linearVelocity;angularVelocity];
end