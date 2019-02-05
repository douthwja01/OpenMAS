function b = bernsteincoeff(p)
%BERNSTEINCOEFF    Bernstein coefficients for polynomial p over [0,1]
%
%   b = bernsteincoeff(p)
%
%For p denoting a polynomial of degree n and B(n,i) denoting the i-th Bernstein 
%polynomial of degree n it is
%
%  p = sum( i=0:n , b(i)*B(n,i)(x) ) .
%
%For convenience, the result b is a polynomial, so that b(i) denotes the i-th Bernstein
%coefficient for i=0:n.
%
%Note that the convex hull of all ( x , p(x) ) for x in [0,1] is contained in the
%convex hull of { ( i/n , b(i) ) : 0<=i<=n } .
%Especially, b(0)=p(0), b(1)=p(1) and the infsup(min(b.c),max(b.c)) contains the range 
%of p over [0,1].
%
%To compute Bernstein coefficients with respect to another interval [a,b] use
%  q = ptrans(p,a,b,0,1) instead of p.
%

% written  10/25/02     S.M. Rump
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

  k = size(p.e,2);          % number of variables
  b = p;
  if k==1                   % univariate case
    
    n = length(p.c)-1;
    for i=0:n
      j = 0:i;
      b.c(n-i+1) = sum(binom(i,j)./binom(n,j).*p.c(n-j+1));
    end
    
  else                      % multivariate case
    
    emax = max(p.e,[],1);   % vector of dimensions
    emax1 = emax+1;
    N = prod(emax1);
    
    % generate exponents
    b.e = (0:emax(1))';
    for i=2:k
      c = ones(prod(emax1(1:i-1)),1)*(0:emax(i));
      b.e = [ repmat(b.e,emax1(i),1) c(:)];
    end
    
    % generate coefficients
    b.c = typeadj(zeros(N,1),typeof(p.c));
    for j=1:size(p.e,1)     % loop through coefficients of p
      jj = p.e(j,:);        % exponents of coefficient p_j
      I = all( ( b.e >= ones(N,1)*jj ) , 2 );   
      ii = b.e(I,:);        % index set of b_i to be updated 
      s = p.c(j);
      for mu=1:k
        s = s .* ( binom(ii(:,mu),jj(mu)) / binom(emax(mu),jj(mu)) );
      end
      b.c(I) = b.c(I) + s;
    end
    
  end
  
  if rndold
    setround(rndold)
  end
