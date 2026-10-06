% BOUND VALUE
function [bounded_v] = boundValue(v,lower,upper)
% This function simply bounds a value between two given values

% Input sanity checks
assert(size(lower,1) == size(upper,1),'Please specify the same number of upper and lower bounds');

% Determine behaviour
switch size(lower,1)
    case size(v,1)      % One bound per dimension
        % Apply each bound 
        for ind = 1:size(v,1)
            if v(ind) < lower(ind)
                v(ind) = lower(ind);
            end
            if v(ind) > upper(ind)
                v(ind) = upper(ind);
            end
        end
    case 1              % A single uni-lateral bound
        % Apply bound to all
        for ind = 1:size(v,1)
            if v(ind) < lower
                v(ind) = lower;
            end
            if v(ind) > upper
                v(ind) = upper;
            end
        end
    otherwise
        error('Input not valid, please provide either a bound per dimension, or a unilateral value.');
end
bounded_v = v;  % Assign new value
end