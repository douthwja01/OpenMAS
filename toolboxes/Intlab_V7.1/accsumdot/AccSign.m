function S = AccSign(p)
%ACCSIGN      Computes sign(sum(p))
%
%   S = AccSign(p)
%
%On return, S is the sign of sum(p), also in the presence
%  of underflow. Input vector p may be single or double precision.
%
%Implements Algorithm 8.2 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008. 
%Requires (4m+2)n flops for m executions of repeat-until loop in Transform.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ~isreal(p)
    error('AccSign for real input only')
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isa(p,'double')
    nmax = 2^52;            % nmax = 450,359,962,737,0496
  else
    nmax = 2^23;            % nmax = 8,388,608
  end
  if length(p)>nmax
    error(['maximum length of input vector for AccSign ' int2str(nmax) '.'])
  end

  kPhi = 1;
  [tau1,tau2,p] = Transform(p,0,kPhi);
  S = sign(tau1);

  if rndold
    setround(rndold)
  end
  