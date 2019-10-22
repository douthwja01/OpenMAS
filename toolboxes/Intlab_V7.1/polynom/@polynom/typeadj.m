function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  08/28/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% polynom\@polynom\typeadj:  a  must be polynom
  switch TYPE
    case 'polynom',        c = mid(a);
    case 'polynomintval',  c = intval(a);
  otherwise
    error('invalid type in call of typeadj')
  end
