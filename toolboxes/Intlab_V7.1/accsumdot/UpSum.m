function res = UpSum(p)
%UPSUM        Rounded upwards result of sum(p)
%
%   resD = DownSum(p)
%
%On return, resD is sum(p) rounded upwards, also in the presence
%  of underflow. Input vector p may be single or double precision.
%
%Adapted from Algorithm 7.5 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%Requires (4m+4k+4)n flops for m,k executions of the repeat-until loops
%  in the calls of TransformK, respectively.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ~isreal(p)
    res = complex(UpSum(real(p)),UpSum(imag(p)));
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
    nmax = 2^26-2;                          % nmax = 67,108,864
  else
    nmax = 2^12-2;                          % nmax = 4,094
  end

  if length(p)<=nmax                        % not huge dimension
    [res,R,p,sigma,Ms] = TransformK(p,0);   % s-res = R+sum(p)
    delta = TransformK(p,R,sigma,Ms);       % delta faithful rounding of s-res
  else                                      % huge dimension
    [res,tau1,tau2,tau,p] = AccSumHugeN(p); % s=tau1+tau2+sum(tau)+sum(p)
    delta = AccSumHugeN([tau1;tau2;tau(:);p(:);-res]);
  end
  if delta>0                                % s > res
    res = succ(res);
  end
  
  if rndold
    setround(rndold)
  end
