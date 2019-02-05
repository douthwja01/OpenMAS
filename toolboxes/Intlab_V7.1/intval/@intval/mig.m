function c = mig(a);
%MIG          mignitude(a)  for interval scalars, vectors and matrices
%
%   c = mig(a)
%

% written  11/09/98     A. Neumaier
% modified 09/02/00     S.M. Rump   rounding unchanged after use
% modified 09/10/02     S.M. Rump   improvement in performance and for complex case
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of huge arrays
% modified 01/06/05     S.M. Rump   semicolon added
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if a.complex
    if isequal(a.rad,0)
      c = abs(a.mid);
    else
      setround(-1)
      c = abs(a.mid)-a.rad;
      setround(0)
    end
    c = c .* ( c>=0 );
  else
    infa = a.inf;
    supa = a.sup;
    c = min(abs(infa),abs(supa)) .* ( ~( ( infa<=0 ) & ( supa>=0 ) ) );
  end
  
  if rndold
    setround(rndold)
  end
