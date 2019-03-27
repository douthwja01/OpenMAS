% THE QUATERNION DIFFERENCE
function [dq] = qDifference(q1,q2)
% This function gets the quaternion that will rotate from the
% q1 orientation to the q2 orientation [q1->q2].

assert(size(q1,1) == 4,'The initial quaternion (q1) must be a 4x1 column vector')
assert(size(q1,1) == 4,'The terminal quaternion (q2) must be a 4x1 column vector')

% Get the quaternion inverse
[q1_inv] = qInverse(q1);
% Calculate difference
dq = qMultiply(q1_inv,q2);
% Re-normalise
dq = unit(dq);
end