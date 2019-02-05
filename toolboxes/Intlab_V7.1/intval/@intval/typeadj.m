function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  10/16/98     S.M. Rump
% modified 12/06/98     S.M. Rump
% modified 12/18/02     S.M. Rump  Hessians added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% intval\@intval\typeadj:  a  must be intval
  switch TYPE
    case 'double',         c = mid(a);
    case 'intval',         c = a;
    case 'gradient',       c = gradient(mid(a));
    case 'gradientintval', c = gradient(a);
    case 'hessian',        c = hessian(mid(a));
    case 'hessianintval',  c = hessian(a);
    case 'taylor',         c = taylor(mid(a));
    case 'taylorintval',   c = taylor(a);
    case 'slope',          c = slope(a);
    case 'polynom',        c = polynom(mid(a));
    case 'polynomintval',  c = polynom(a);
  otherwise
    error('invalid type in call of typeadj')
  end
