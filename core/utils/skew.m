% SKEW MATRIX
function [V] = skew(v)
assert(numel(v) == 3,'Input vector must be [3x1]');
% SKEW THE VECTOR COMPONENTS
V = [ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end