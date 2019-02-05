function i = end(x,k,n)
%END          Overloaded functions end, specifies last index
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if n==1           % call as one-dimensional array
    i = length(x.x(:));
  else
    i = size(x.x,k);
  end
