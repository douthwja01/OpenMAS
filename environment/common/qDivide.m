% QUATERNION DIVISION
function [qDiv] = qDivide(q,r)
% The function compute the division of quaternion q, by
% quaterion v.
% Associated block:
% "Quaternion Division"

qDiv = zeros(4,1);
% Normalise the second vector
qDiv(1) =  q(1)*r(1) + q(2)*r(2) + q(3)*r(3) + q(4)*r(4);
qDiv(2) = -q(1)*r(2) + q(2)*r(1) + q(3)*r(4) - q(4)*r(3);
qDiv(3) =  q(1)*r(3) + q(2)*r(4) + q(3)*r(1) + q(4)*r(2);
qDiv(4) = -q(1)*r(4) + q(2)*r(3) - q(3)*r(2) + q(4)*r(1);
qDiv = qDiv./(norm(r));
end