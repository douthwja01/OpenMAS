% NUMERICALLY GET THE ROTATIONS BETWEEN TWO TRIADS
function [q] = getNumericalTriadRotation(referenceTriad,targetTriad)
% This function defines the quaternion nescessary to rotate one
% axis triad into alignment with a second triad. This is done by
% solving the rotation matrix-defined symbolic expressions for
% the equivalent quaternion rotation.

syms q1 q2 q3 q4 real

% Define the symbolic quaternion rotation matrix
R_q = zeros(3);
R_q(1,1) = q1^2 + q2^2 - q3^2 - q4^2;
R_q(1,2) = 2*(q1*q4 + q2*q3);
R_q(1,3) = 2*(q2*q4 - q1*q3);
R_q(2,1) = 2*(q3*q2 - q1*q4);
R_q(2,2) = q1^2 - q2^2 + q3^2 - q4^2;
R_q(2,3) = 2*(q1*q2 + q3*q4);
R_q(3,1) = 2*(q1*q3 + q2*q4);
R_q(3,2) = 2*(q3*q4 - q1*q2);
R_q(3,3) = q1^2 - q2^2 - q3^2 + q4^2;
% DEFINE THE SYMBOLIC EXPRESSIONS
expressions = [R_q*referenceTriad(:,1) - targetTriad(:,1);
    R_q*referenceTriad(:,2) - targetTriad(:,2);
    R_q*referenceTriad(:,3) - targetTriad(:,3);
    q1^2 + q2^2 + q3^2 + q4^2 - 1];
expressions = symfun(expressions,[q1 q2 q3 q4]);                % Define symbolic function in q_{1-4}
modeqn = @(x)double(expressions(x(1),x(2),x(3),x(4)));          % Rearrange to be in terms of q1,q2 etc
% DEFINE SOLVER OPTIONS
OPTIONS = optimoptions('fsolve','Algorithm','levenberg-marquardt');
% PASS SYMBOLIC EQUATIONS TO F-SOLVE
x0 = [ 1 0 0 0];                 % INITIAL QUATERNION GUESS
q = fsolve(modeqn,x0,OPTIONS);   % SOLVE FOR THE ROTATION QUATERNION
% RE-NORMALISE THE QUATERNION
[q] = OMAS_geometry.qUnit(q);
end