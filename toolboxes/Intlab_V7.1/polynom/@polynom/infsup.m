function infsup(p)
%INFSUP       Display of polynom in "infsup" mode
%

% written  07/22/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  olddisp = intvalinit('display',0);
  intvalinit('displayinfsup',0);
  display(p,inputname(1));
  intvalinit(olddisp,0);
