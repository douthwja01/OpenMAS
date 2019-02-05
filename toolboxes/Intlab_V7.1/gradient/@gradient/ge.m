function res = ge(a,b)
%GE           Implements  a >= b  for gradients, compares only a.x and b.x
%

% written  10/16/98     S.M. Rump
% modified 12/18/02     S.M. Rump  complex comparison 
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isa(a,'gradient')
    if isreal(a) & isreal(b.x)
      res = ( a>=b.x );
    else
      res = ( real(a)>=real(b.x) ) & ( imag(a)>=imag(b.x) );
    end
  elseif ~isa(b,'gradient')
    if isreal(a.x) & isreal(b)
      res = ( a.x>=b );
    else
      res = ( real(a.x)>=real(b) ) & ( imag(a.x)>=imag(b) );
    end
  else
    if isreal(a.x) & isreal(b.x)
      res = ( a.x>=b.x );
    else
      res = ( real(a.x)>=real(b.x) ) & ( imag(a.x)>=imag(b.x) );
    end
  end
