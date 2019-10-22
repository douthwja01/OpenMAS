function Ac = compmat(A);
%COMPMAT      Comparison matrix for interval matrices
%
%   Ac = compmat(A)
%

% written  11/09/98     A. Neumaier
% modified 08/07/02     S.M. Rump    abss instead of abs
% modified 04/04/04     S.M. Rump    set round to nearest for safety
% modified 04/06/05     S.M. Rump    rounding unchanged
% modified 09/28/08     S.M. Rump    check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  Ac = -mag(A);
  Ac = ( Ac - diag(diag(Ac)) ) + diag(mig(diag(A))) ;
  
  if rndold
    setround(rndold)
  end
