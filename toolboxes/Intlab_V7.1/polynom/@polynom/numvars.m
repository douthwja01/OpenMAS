function k = numvars(p)
%NUMVARS      Number of variables of polynomial
%
%   k = numvars(p)
%
%Careful: number of variables of a constant polynomial is zero.
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if isempty(p.v)
    k = 0;
  else
    k = size(p.e,2);
  end
