function a = mpower(a,b)
%MPOWER       Implements  a ^ b  for slopes (b must be scalar)
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if prod(size(b))~=1
    error('slope mpower only for scalar exponent')
  end
  a = a .^ b;
