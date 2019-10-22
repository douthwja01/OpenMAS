function R = InvIllco(A,iter)
%INVILLCO     Inverse of extremely ill-conditioned matrices
%
%   R = InvIllco(A,iter)
%
%On return, R is a cell array such that sum(R{i}) is an approximate inverse
%  of A and  I-R*A  is convergent. 
%The input matrix A itself may be a cell array. This is convenient to store
%  not exactly representable input data in higher precision. Then R is an
%  approximate inverse of sum(A{i}).
%
%The parameter iter is optional. If specified, an extra iteration is 
%  executed to produce  I-R*A  of order eps.
%
%Example:
%
%  n = 20;              % dimension os matrix
%  C = 1e50;            % anticipated condition number
%  A = randmat(n,C);    % ill-conditioned matrix
%  R = InvIllco(A);     % approximate inverse, stored in cell array R
%  norm(AccDot(R,A,-1,eye(n)),'fro')     
%
%Implements algorithm InvIllco from
%  S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
%    Japan J. Indust. Appl. Math. (JJIAM), 26:249-277, 2009.
%
%Reference implementation! Slow due to interpretation!
%

% written  06/23/08     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, warning
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  wng = warning;
  warning off

  if iscell(A)
    n = size(A{1},1);
    R = inv(A{1});
  else
    n = size(A,1);
    R = inv(A);
  end
  k = 1;
  iter = ( nargin==2 );
  while 1
    k = k+1;
    P = ProdKL(R,A,k,1);
    X = inv(P);
    while any(any(isinf(X)))
      X = inv(P.*(1+eps*randn(n)));
    end
    R = ProdKL(X,R,k,k);
    if norm(X,'fro')*norm(P,'fro')<.01/eps
      if iter, iter = 0; else break, end
    end
  end
  
  if rndold
    setround(rndold)
  end
  warning(wng)
