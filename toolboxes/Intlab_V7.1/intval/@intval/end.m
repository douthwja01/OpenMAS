function i = end(x,k,n)
%END          Overloaded functions end, specifies last index
%

% written  10/03/02     S.M. Rump
%

  if x.complex
    xx = x.mid;
  else
    xx = x.inf;  
  end
  if n==1           % call as one-dimensional array
    i = length(xx(:));
  else
    i = size(xx,k);
  end
