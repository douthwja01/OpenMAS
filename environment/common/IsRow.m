
function [flag] = IsRow(v,n)
% Default setting
flag = false;

% Get vector dimensions
[a,b] = size(v);

% Check it is a column
if a == 1 && nargin == 1
    flag = true;
    return;
end  

% Otherwise check its length against n
if b == n           % Is not same length
    flag = true;
    return       
end
end
