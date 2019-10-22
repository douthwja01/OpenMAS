function [res,tau1,tau2,tau,p] = AccSumHugeN(p)
%ACCSUMHUGEN  Faithful rounding of sum(p) for huge dimension
%
%   res = AccSumHugeN(p)
%
%On return, res is a faithful rounding of sum(p), also in the presence
%  of underflow. Input vector p may be single or double precision.
%
%Implements Algorithm 8.1 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%Requires (4m+4K+3)n flops for m executions of repeat-until loop in 
%  Transform and K executions of the while-loop.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ~isreal(p)
    res = complex(AccSumHugeN(real(p)),AccSumHugeN(imag(p)));
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
    nmax = 2^50;            % nmax = 1,125,899,906,842,624
  else
    nmax = 2^21;            % nmax = 2,097,152
  end
  if length(p)>nmax
    error(['maximum length of input vector for AccSumHugeN ' int2str(nmax) '.'])
  end
 
  kPhi = 2;
  [tau1,tau2,p,sigma,Ms] = Transform(p,0,kPhi);   % Ms = 2^M
  if sigma<=realmin             % p_i identical zero
    res = tau1;
    tau = 0*res;                % make sure tau has same precision
    if rndold, setround(rndold); end
    return
  end
  if isa(p,'double'), prec='double'; else prec='single'; end    
  u = 0.5*eps(prec);
  tau = zeros(1,16);
  K = 0;
  phi = Ms*u;
  factor = 2*Ms*Ms*u;
  sigmas = phi*sigma;
  while 1
    K = K+1;
    sigma = sigmas;
    [tau(K),p] = ExtractVector(p,sigma);
    sigmas = phi*sigma;
    if ( factor*sigma<=abs(tau1) ) | ( sigma<=realmin )
      taus = tau2 + sum(p);
      for k=K:-1:1
        taus = taus + tau(k);
      end
      res = tau1 + taus;
      if rndold, setround(rndold); end
      return
    end
  end
  