function r = isnan(a)
%ISNAN        Array of 1's for NaN components
%
%   r = isnan(a)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/17/06     S.M. Rump  huge arrays and a.rad=0
%

  if a.complex
    if isequal(a.rad,0)
      r = isnan(a.mid);
    else
      r = isnan(a.mid) | isnan(a.rad) ;
    end
  else
    r = isnan(a.inf) | isnan(a.sup) ;
  end
