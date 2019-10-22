function c = typeadj(a,TYPE)
%TYPEADJ      typecast of  a  to type TYPE
%
%   c = typeadj(a,TYPE)
%
%Adjust type of  a  to be TYPE. For possible values of TYPE, see intval\typeof.
%Typical application:
%   c = typeadj(a,typeof(x));
%then  c = a  but with c being of type of x.
%Adjustment of c is as follows:
%
%    typeof(a)\TYPE  |   d     i      g      gi       h      hi      s      p       pi     t     ti
%  --------------------------------------------------------------------------------------------------
%         d          |   *   i(a)    g(a)  g(i(a))   h(a)  h(i(a))  s(a)   p(a)   p(i(a)) t(a)  t(i(a))
%         i          |  m(a)   *   g(m(a))   g(a)  h(m(a))  h(a)    s(a)  p(m(a))  p(a)  t(m(a)) t(a) 
%         g          |   -     -      *      i(a)     -      -       -      -       -      -     -
%         gi         |   -     -     m(a)     *       -      -       -      -       -      -     -
%         h          |   -     -      -       -      h(a)  h(i(a))   -      -       -      -     -
%         hi         |   -     -      -       -      m(a)    *       -      -       -      -     -
%         s          |   -     -      -       -       -      -       *      -       -      -     -
%         p          |   -     -      -       -       -      -       -      *    p(i(a))   -     -
%         pi         |   -     -      -       -       -      -       -     m(a)    p(a)    -     -
%         t          |   -     -      -       -       -      -       -      -       -      *     -
%         ti         |   -     -      -       -       -      -       -      -       -     m(a)   *
%
%with   *    c = a;
%      m(a)  mid(a)
%      i(a)  intval(a)
%      g(a)  gradient(a)
%      h(a)  hessian(a)
%      s(a)  slope(a)
%      p(a)  polynom(a)
%      t(a)  taylor(a)
%       -    not applicable (error)
%

% written  10/16/98     S.M. Rump
% modified 12/06/98     S.M. Rump
% modified 12/18/02     S.M. Rump  Hessians added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/22/09     S.M. Rump  taylor added
%

% intval\typeadj:  a  must be double
  switch TYPE
    case 'double',         c = a;
    case 'intval',         c = intval(a);
    case 'gradient',       c = gradient(a);
    case 'gradientintval', c = gradient(intval(a));
    case 'hessian',        c = hessian(a);
    case 'hessianintval',  c = hessian(intval(a));
    case 'slope',          c = slope(a);
    case 'polynom',        c = polynom(a);
    case 'polynomintval',  c = polynom(intval(a));
    case 'taylor',         c = taylor(a);
    case 'taylorintval',   c = taylor(intval(a));
  otherwise
    error('invalid type in call of typeadj')
  end
  