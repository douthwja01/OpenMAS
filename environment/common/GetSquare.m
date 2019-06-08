function [coords] = GetSquare(x1,x2,x3,x4)

% This function draws a square defined by vertices

% Sanity check
assert(numel(x1) <= 3,"Coordinates can have upto three dimensions.");
assert(iscolumn(x1) == iscolumn(x2) == iscolumn(x3) == iscolumn(x3),"The coordinate vectors must be the same length.");

coords = zeros(5,3);
for i = 1:numel(x1)
    coords(:,i) = [x1(i);x2(i);x3(i);x4(i);x1(i)];
end

end