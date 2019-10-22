
function [n_ccw] = CCWnormal_2D(v)
    % Input sanity check
    assert(size(v,1) == 2,"Expecting a 2D column vector.");
    % Simply get the CCW normal
    n_ccw = [-v(2);v(1)];
end

% C++ implementation
% public static Vector2 PerpendicularCounterClockwise(this Vector2 vector2)
% {
%     return new Vector2(-vector2.Y, vector2.X);
% }