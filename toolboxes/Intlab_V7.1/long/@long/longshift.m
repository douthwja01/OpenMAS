function A = longshift(A,r)
%LONGSHIFT    Shifts long number A by r bits
%
%   A = longshift(A,r)
%
%The same as A * long(2)^r
%

% written  01/10/04     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_LONG_LOGBETA = getappdata(0,'INTLAB_LONG_LOGBETA');

  q = floor(r/INTLAB_LONG_LOGBETA);
  A = struct(A);
  A.exponent = A.exponent + q;
  rem = r - q*INTLAB_LONG_LOGBETA;
  A.mantissa = pow2(A.mantissa,rem);
  A = normalize(A);
  A = normalizefirst(A);
  A = class(A,'long');
    
  if rndold
    setround(rndold)
  end
