% GET THE QUATERNION DIFFERNTIAL (as a function of omega)
function [q_dot] = qDifferential(q0,omega)
% This function was created to allow the quaternion
% differential to be called seperately to in the integration
% function. This allows the quaternion 'q0' to be integrated
% externally.

% Calculate the integration drift normalising elements
K = 1; % k*dt <= 1
lambda = 1 - (q0(1)^2 + q0(2)^2 + q0(3)^2 + q0(4)^2);

% Calculate the quaternion difference                          % NOTATION B
q_dot = 0.5*[ K*lambda, -omega(1), -omega(2), -omega(3);
              omega(1),  K*lambda,  omega(3), -omega(2);
              omega(2), -omega(3),  K*lambda,  omega(1);
              omega(3),  omega(2), -omega(1),  K*lambda]*q0;
end