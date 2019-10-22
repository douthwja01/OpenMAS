
function [n_cw] = CWnormal_2D(v)
    % Input sanity check
    assert(size(v,1) == 2,"Expecting a 2D column vector.");
    % Simply get the C normal
    n_cw = [v(2);-v(1)];
end


% C++ implementation
% public static Vector2 PerpendicularClockwise(this Vector2 vector2)
% {
%     return new Vector2(vector2.Y, -vector2.X);
% }