function r = isinf(a)
%ISINF        Array of 1's for inf components
%
%   r = isinf(a)
%

% written  08/07/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/17/06     S.M. Rump  huge arrays and a.rad=0
%

  if a.complex
    if isequal(a.rad,0)
      r = isinf(a.mid);
    else
      r = isinf(a.mid) | isinf(a.rad) ;
    end
  else
    r = isinf(a.inf) | isinf(a.sup) ;
  end
