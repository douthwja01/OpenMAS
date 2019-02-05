function res = SumKL(p,K,L)
%SUMKL        Summation 'as if' computed in K-fold precision and stored in L results
%
%   res = SumKL(p,K,L)
%
%On return, sum(res) approximates sum(p) with accuracy as if computed 
%  in K-fold precision, where res comprises of L elements. 
%  Default for L is 1.
%
%Implements algorithm SumKL from
%  S.M. Rump: Inversion of extremely ill-conditioned matrices in floating-point,
%      Japan J. Indust. Appl. Math. (JJIAM), 26:249-277, 2009.
%
%Reference implementation! Slow due to interpretation!
%

% written  06/23/08     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin==2
    L = 1;
  end
  
  n = length(p);
  for i=1:K-L
    p = VecSum(p);
  end
  res = zeros(1,L);
  for k=0:L-2
    p(1:n-k) = VecSum(p(1:n-k));
    res(k+1) = p(n-k);
  end
  res(L) = sum(p(1:n-L+1));
  
  if rndold
    setround(rndold)
  end
