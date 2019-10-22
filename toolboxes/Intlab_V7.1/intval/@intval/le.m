function res = le(a,b)
%LE           Implements  a <= b  elementwise for intervals a and b
%
%  if true,  a  is definitely less than or equal to  b
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isa(a,'intval')
    a = intval(a);
  end
  if ~isa(b,'intval')
    b = intval(b);
  end

  if a.complex | b.complex
    res = real(sup(a)) <= real(inf(b)) & imag(sup(a)) <= imag(inf(b)) ;
  else
    res = sup(a) <= inf(b) ;
  end
