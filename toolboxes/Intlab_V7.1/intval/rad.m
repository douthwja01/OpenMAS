function r = rad(a)
%RAD          Implements  rad(a)
%
%   r = rad(a)
%

%Implementation for point data
%

% written  10/16/98     S.M. Rump
% modified 06/01/99     S.M. Rump  zeros of proper size
% modified 10/13/99     S.M. Rump  rad(NaN) = NaN
% modified 04/03/04     S.M. Rump  sparse case
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/17/12     S.M. Rump  treatment of NaN
%

  if issparse(a)
    [m n] = size(a);
    r = sparse([],[],[],m,n,0);
  else
    r = zeros(size(a));
  end
  
  if issparse(a)    % take care of huge arrays
    [i,j] = find(isnan(a));
    [m n] = size(a);
    r = r + sparse(i,j,NaN,m,n);
  else
    r(isnan(a)) = NaN;
  end
