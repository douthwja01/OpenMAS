% GET THE NORM OF AN INTERVAL VECTOR
function [intNorm] = intNorm(intVector,zeroOmitted)
% INPUT HANDLING
if ~exist('zeroOmitted','var')
    zeroOmitted = 0;
end

if ~isintval(intVector)
    % If the input value is not an interval
    intNorm = norm(intVector);
    return
end

% IF ZERO NOT IN RANGE
if zeroOmitted
    for i = 1:size(intVector)
        s = 0; % Collect the square terms
        for j = 1:size(intVector)
            if i == j
                continue
            end
            s = s + intVector(j)^2;
        end
        x = intVector(i);
        
        % COMPUTE THE INTERVAL NORM WHERE 0 IS NOT INCLUDED
        n(i) = sign(inf(x))/sqrt(1+(s/x^2));
    end
else
    % ZERO IS PART OF THE RANGE ( GENERALISED )
    for i = 1:size(intVector,1)
        s = 0; % Collect the square terms
        for j = 1:size(intVector,1)
            if i == j
                continue
            end
            s = s + sqr(intVector(j));
        end
        x = intVector(i);
        
        offset = 1E-10;
        upperbound = infsup(offset,inf);
        xi_pos = intersect(x,upperbound);
        lowerbound = infsup(-inf,-offset);
        xi_neg = intersect(x,lowerbound);
        
        ni_upper = 1./sqrt(1 + (s./sqr(xi_pos)));
        ni_lower = -1./sqrt(1 + (s./sqr(xi_neg)));
        
        if isnan(ni_upper)
            n(i) = ni_lower;
        elseif isnan(ni_lower)
            n(i) = ni_upper;
        else
            n(i) = hull(ni_upper,ni_lower);
        end
    end
end
intNorm = n;
% CONFIRM FORMATTING OF THE VECTOR
intNorm = reshape(intNorm,size(intVector));
end