function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% hessian\@hessian\typeadj:  input  a  must be hessian (superior to intval)
  switch TYPE
    case 'hessian',       c = mid(a);
    case 'hessianintval', c = intval(a);
  otherwise
    error('invalid type in call of typeadj')
  end
