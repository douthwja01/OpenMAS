function r = isintval(a)
%ISINTVAL     Returns 1 if  a  is intval hessian
%
%   r = isintval(a)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  r = isa(a.x,'intval');
