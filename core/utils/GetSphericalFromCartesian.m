% CONVERT CARTESIAN TO SPHERICAL
function [r,phi,theta]  = GetSphericalFromCartesian(p)
% INPUTS:
% p - Cartesian point [x,y,z]
% OUTPUTS:
% r   - Range (m)
% phi - XY angle (rad)
% theta - XZ angle (rad)

r = norm(p);                                        % The range
phi = atan2(p(2),p(1));                             % The angle made in the azimuth (bearing)
theta = atan2(p(3),sqrt(p(1).^2 + p(2).^2));        % The elevation (vertical bearing)
end