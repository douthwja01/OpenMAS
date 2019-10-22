function r = minus(p,q)
%MINUS        Polynomial subtraction  p - q
%
%p or q may be scalar (interval)
%

% written  09/02/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = p + (-q) ;
