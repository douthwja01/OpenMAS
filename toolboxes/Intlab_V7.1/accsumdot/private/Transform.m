function [tau1,tau2,p,sigma,Ms] = Transform(p,rho,kPhi,sigma,Ms)
%TRANSFORM    Error-free transformation of rho+sum(p)
%
%   [tau1,tau2,ps] = Transform(p,rho,kPhi)
%
%On return, tau1+tau2+sum(ps_i) = rho+sum(p_i), also in the presence
%  of underflow. Input vector p may single or double precision.
%
%Parameter rho optional, default 0
%
%Parameter kPhi (default 0) specifies stopping criterion:
%  kPhi  0   2^(2M) eps     for AccSum, TransformK, AccSumK, NearSum
%        1   2^M eps        for AccSign
%        2   eps            for AccSumHugeN
%Length of input vector not checked.
%
%Implements Algorithm 4.1 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, SIAM J. Sci. Comput., 31(1):189-224, 2008.
%and Algorithm 3.3 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%Requires (4m+2)n flops for m executions of the repeat-until loop if
%  sigma,Ms not given, otherwise 4mn flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  if isa(p,'double'), prec='double'; else prec='single'; end
  u = 0.5*eps(prec);
  if nargin<5                           % standard call
    if nargin<2, rho=0; end
    if nargin<3, kPhi=0; end
    n = length(p);                      % initialization
    mu = max(abs(p));                   % abs(p_i) <= mu
    if ( n==0 ) | ( mu==0 )             % no or only zero summands
      tau1 = rho;
      tau2 = 0;
      sigma = 0;
      Ms = 0;
      return
    end
    Ms = NextPowerTwo(n+2);             % n+2 <= 2^M = Ms
    sigma = Ms*NextPowerTwo(mu);        % first extraction unit
    if ~isfinite(sigma)
      error('overflow occurred in Transform')
    end
  else                                  % sigma and Ms already known
    sigma = Ms*u*sigma;
  end    
  phi = Ms*u;                           % factor to decrease sigma    
  switch kPhi                           % choose stopping criterion
    case 0, Phi = Ms*Ms*u;              % AccSum, TransformK, AccSumK, NearSum
    case 1, Phi = Ms*u;                 % AccSign
    case 2, Phi = 8*Ms*u;               % AccSumHugeN
  end
  t = rho;
  while 1
    q = ( sigma + p ) - sigma;          % [tau,p] = ExtractVector(sigma,p);
    tau = sum(q);                       % sum of leading terms
    p = p - q;                          % remaining terms
    tau1 = t + tau;                     % new approximation
    if ( abs(tau1)>=Phi*sigma ) | ( sigma<=realmin(prec) )      
      tau2 = tau - ( tau1 - t );        % [tau1,tau2] = FastTwoSum(t,tau)
      return
    end
    t = tau1;                           % sum t+tau exact
    if t==0                             % accelerate case sum(p)=0
      [tau1,tau2,p,sigma,Ms] = Transform(p,0,kPhi);  % sum of remainder part
      return
    end
    sigma = phi*sigma;                  % new extraction unit
  end
  