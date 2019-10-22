function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  05/22/09     S.M. Rump
%

% taylor\@taylor\typeadj:  a  must be taylor (superior to intval)
  switch TYPE
    case 'taylor',       c = mid(a);
    case 'taylorintval', c = a;
  otherwise
    error('invalid type in call of typeadj')
  end
