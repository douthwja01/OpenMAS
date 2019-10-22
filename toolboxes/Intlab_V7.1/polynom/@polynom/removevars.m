function p = removevars(p)
%REMOVEVARS   Remove superfluous variables
%
%   q = removevars(p)
%
%Variables with only zero coefficients are removed. If result is
%  constant, set q.v of variables is empty.
%

% written  09/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ischar(p.v)                  % univariate polynomial
    if p.e==0                     % constant polynomial
      p.v = '';
    end
  else
    index = ( sum(p.e,1) ~= 0 );  % index set of actually needed variables
    if ~any(index)                % constant polynomial
      p.e = 0;
      p.v = '';
    elseif ~all(index)
      % multivariate, but some variables superfluous
      p.e = p.e(:,index);
      p.v = p.v(index);
      if size(p.e,2)<=1           % result univariate
        n = max(p.e);             % p.e is vector
        c = p.c;
        p.c = typeadj( zeros(1,n+1) , typeof(c) );
        p.c(n+1-p.e) = c;
        p.e = n;
        p.v = p.v{1};             % result cannot be constant (covered above)
      end
    end
  end
  
  if rndold
    setround(rndold)
  end
