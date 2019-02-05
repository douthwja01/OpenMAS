function A = bezout(p,q)
%BEZOUT       Bezout matrix of two univariate polynomials p and q
%
%   A = bezout(p,q)
%
%Both polynomials must be univariate with the same independent variable.
%The (symmetric) Bezout matrix is singular iff p and q have a root in common.
%For p or q being interval polynomials, A is an interval matrix.
%The call
%
%   A = bezout(p)
%
%is the same as   A = bezout(p,p') .
%The definition follows Fiedler, Special matrices and their applications.
%

% written  02/03/04     S.M. Rump
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

  if nargin==1
    q = p';
  end

  if ( size(p.v)>1 ) | ( size(q.v)>1 )
    error('both polynomials must be univariate')
  end
  if p.v~=q.v
    error('both polynomials must depend on the same variable')
  end
  
  np = p.e;
  nq = q.e;
  n = max(np,nq);

  A = toeplitz([p.c(np+1) zeros(1,n-1)],[fliplr(p.c) zeros(1,2*n-np-1)]);
  B = toeplitz([q.c(nq+1) zeros(1,n-1)],[fliplr(q.c) zeros(1,2*n-nq-1)]);
  A = flipud( A(:,n+1:end)*B(:,1:n) - B(:,n+1:end)*A(:,1:n) );
  
  if rndold
    setround(rndold)
  end
