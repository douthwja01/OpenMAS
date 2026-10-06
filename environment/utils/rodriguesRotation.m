% ROTATE A VECTOR ABOUT AN AXIS VECTOR (RODRIGUES)
function [v_rotated] = rodriguesRotation(u,k,theta)
% v - Vector to be rotated
% k - Is the rotation axis
% Theta - The angle the vector is to be rotated through

assert(numel(u) == 3,'Rotation vector must be of size [3x1].')
assert(numel(k) == 3,'The rotation axis must be of size [3x1]');
assert(numel(theta) == 1,'The rotation angle %.0f must be a scalar',theta);

[m,n] = size(u);
if (m ~= 3 && n ~= 3)
    error('input vector is/are not three dimensional')
end
if (size(u) ~= size(k))
    error('rotation vector v and axis k have different dimensions')
end

k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
No = numel(u)/3; % number of vectors in array
v_rotated = u; % initialize rotated vector array
if ( n == 3 )
    crosskv = u(1,:); % initialize cross product k and v with right dim.
    for i = 1:No
        crosskv(1) = k(2)*u(i,3) - k(3)*u(i,2);
        crosskv(2) = k(3)*u(i,1) - k(1)*u(i,3);
        crosskv(3) = k(1)*u(i,2) - k(2)*u(i,1);
        v_rotated(i,:) = cos(theta)*u(i,:) + (crosskv)*sin(theta)...
            + k*(dot(k,u(i,:)))*(1 - cos(theta));
    end
else % if m == 3 && n ~= 3
    crosskv = u(:,1); % initialize cross product k and v with right dim.
    for i = 1:No
        crosskv(1) = k(2)*u(3,i) - k(3)*u(2,i);
        crosskv(2) = k(3)*u(1,i) - k(1)*u(3,i);
        crosskv(3) = k(1)*u(2,i) - k(2)*u(1,i);
        v_rotated(:,i) = cos(theta)*u(:,i) + (crosskv)*sin(theta)...
            + k*(dot(k,u(:,i)))*(1 - cos(theta));
    end
end
end