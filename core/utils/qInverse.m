% QUATERNION INVERSE
function [qInv] = qInverse(q)
% The quaternion norm
%             qSqr = q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2;
qSqr = norm(q);
% CALCULATE THE INVERSE OF A GIVEN QUATERNION
qInv = zeros(4,1);
qInv(1) =  q(1)/qSqr;
qInv(2) = -q(2)/qSqr;
qInv(3) = -q(3)/qSqr;
qInv(4) = -q(4)/qSqr;
end