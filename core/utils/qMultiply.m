% QUATERNION MULTIPLICATION
function [qv] = qMultiply(q,v)
% Calculate the product of two quaternions
% Associated block:
% "Quaternion Multiplication"
% Multiply the quaternion elements

assert(size(q,1) == 4 && size(v,1) == 4,...
    'Both quaternion must be provided as 4x1 column vectors')
% Quaternion projection matrix
qv = [v(1), -v(2), -v(3), -v(4);
      v(2),  v(1), -v(4),  v(3);
      v(3),  v(4),  v(1), -v(2);
      v(4), -v(3),  v(2),  v(1)]*q; % Confirmed with matlab website
% Re-normalise the output
qv = unit(qv);
end