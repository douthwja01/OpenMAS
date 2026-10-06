% QUATERNION CONJUGATE
function [qCon] = qConjugate(q)
% This function defines the conjugate of the given quaternion.
% Associated block:
% "Quaternion Conjugate"
qCon = zeros(4,1);
qCon(2) = -q(2);
qCon(3) = -q(3);
qCon(4) = -q(4);
end