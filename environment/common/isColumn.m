
function [flag] = isColumn(v,n)
% Sanity check #1
assert(isnumeric(v),'Expecting a numeric vector.');

% Default setting
flag = false;

% Get vector dimensions
[a,b] = size(v);

% Check it is a column
if b ~= 1
    return
end

% Sanity check #2
if nargin > 1
    assert(isnumeric(n),'Expecting a numeric column length.');
else
    flag = true;
    return
end
% Check its length against n
if a ~= n
    return       % Is not same length
else
    flag = true; % Is same length
end
end