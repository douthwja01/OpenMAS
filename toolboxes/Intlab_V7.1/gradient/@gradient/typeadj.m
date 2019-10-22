function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  10/16/98     S.M. Rump
% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% gradient\@gradient\typeadj:  a  must be gradient (superior to intval)
  switch TYPE
    case 'gradient',       c = mid(a);
    case 'gradientintval', c = a;
  otherwise
    error('invalid type in call of typeadj')
  end
