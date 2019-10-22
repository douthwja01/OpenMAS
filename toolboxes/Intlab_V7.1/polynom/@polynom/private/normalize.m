function r = normalize(r)
%NORMALIZE    normalization of polynomials
%
%For given multivariate polynomial r,
%  On output, zero coefficients are eliminated and identical sets of exponents
%  summed up
%If input polynomial is a constant, then output is constant polynomial.
%

% written  08/31/00     S.M. Rump
%

  % omit zero coefficients
  if size(r.e,2)>1                % multivariate case
    [r.e,r.c] = collect(r.e,r.c);
    index = ( r.c ~= 0 );
    if ~any(index)                % zero polynomial
      r.e = zeros(1,size(r.e,2));
      r.c = typeadj( 0 ,typeof(r.c) );
      return
    elseif ~all(index)
      r.e = r.e(index,:);
      r.c = r.c(index);
    end
  else                            % univariate case
    index = ( r.c==0 );
    if all(index)                 % zero polynomial
      r.e = zeros(1,size(r.e,2));
      r.c = typeadj( 0 , typeof(r.c) );
      return
    end
    [m,i] = min(index);           % first nonzero coefficient in r.c(i)
    if i~=1                       % omit leading zeros
      r.e = r.e-i+1;
      r.c(1:i-1) = [];
    end
  end
