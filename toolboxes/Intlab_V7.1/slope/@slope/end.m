function i = end(x,k,n)
%END          Overloaded functions end, specifies last index
%

% written  10/03/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if n==1           % call as one-dimensional array
    i = size(x.r,1);
  else
    i = x.size(k);
  end
