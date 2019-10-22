function i = end(x,k,n)
%END          Overloaded functions end, specifies last index
%

% written  10/03/02     S.M. Rump
%

  if n==1           % call as one-dimensional array
    i = length(x.x(:));
  else
    i = size(x.x,k);
  end
