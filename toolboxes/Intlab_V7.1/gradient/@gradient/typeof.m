function t = typeof(a)
%TYPEOF       Type of a
%
%   t = typeof(a)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  10/16/98     S.M. Rump
% modified 12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% gradient\@gradient\typeof:  a  must be gradient
  if isa(a.x,'intval')
    t = 'gradientintval';
  else
    t = 'gradient';
  end
