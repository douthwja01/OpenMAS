function A = sylvester(p,q)
%SYLVESTER    Sylvester matrix of two univariate polynomials p and q
%
%   A = sylvester(p,q)
%
%Both polynomials must be univariate with the same independent variable.
%The Sylvester matrix is singular iff p and q have a root in common.
%For p or q being interval polynomials, A is an interval matrix.
%The call
%
%   A = sylvester(p)
%
%is the same as   A = sylvester(p,p') .
%

% written  08/09/02     S.M. Rump
% modified 02/03/04     S.M. Rump  one input allowed
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

  if np==0
    if nq==0
      A = p.c*q.c;
    else
      A = p.c*eye(nq);
    end
  elseif nq==0
    A = q.c*eye(np);
  else
    A = [ toeplitz([p.c(1) zeros(1,nq-1)],[p.c zeros(1,nq-1)]) ; ...
          toeplitz([q.c(1) zeros(1,np-1)],[q.c zeros(1,np-1)]) ];
  end
  
  if rndold
    setround(rndold)
  end
