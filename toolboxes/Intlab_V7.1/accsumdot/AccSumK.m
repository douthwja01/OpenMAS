function Res = AccSumK(p,K)
%ACCSUMK      Res represents K-fold faithful rounding of sum(p)
%
%   Res = AccSumK(p,K)
%
%On return, Res represents a K-fold faithful rounding of sum(p), also in the 
%  presence of underflow. Input vector p may be single or double precision.
%
%Implements Algorithm 6.4 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%Requires (4m+5K+3)n flops for m executions of repeat-until loop in the
%  first and one execution of the repeat-until loops in subsequent calls
%  of TransformK.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ~isreal(p)
    Res = complex(AccSumK(real(p),K),AccSumK(imag(p),K));
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isa(p,'double')
    nmax = 2^26-2;          % nmax = 67,108,864
  else
    nmax = 2^12-2;          % nmax = 4,094
  end
  if length(p)>nmax
    error(['maximum length of input vector for AccSumK ' int2str(nmax) '.'])
  end

  R = 0;
  if isa(p,'double'), prec='double'; else prec='single'; end
  Res = zeros(1,K,prec);
  [Res(1),R,p,sigma,Ms] = TransformK(p,R);
  if abs(Res(1))<=realmin(prec)
    if rndold, setround(rndold); end
    return
  end
  for k=2:K
    [Res(k),R,p,sigma,Ms] = TransformK(p,R,sigma,Ms);
    if abs(Res(k))<=realmin(prec)
      if rndold, setround(rndold); end
      return
    end
  end
  
  if rndold
    setround(rndold)
  end
