function [s,e,m] = splitdble(d)
%SPLTDBLE     Split double vector into sign, mantissa and exponent
%
%  [s,e,m] = splitdble(d)
%
%such that  d = s * m * 2^e  for integers s,m,e
%fixes gradual underflow bug
%

% written  12/30/98     S.M. Rump
%

  s = sign(d);
  s(s==0) = 1;

  d = abs(d);
  [m,e] = log2(d);

  index = ( d<realmin );
  if any(index)
    [m1,e1] = log2(d(index)*2^50);
    m(index) = m1;
    e(index) = e1 - 50;
  end
