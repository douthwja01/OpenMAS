function res = PrecSum(p,K)
%PRECSUM      Fast summation better than K-fold precision
%
%   res = PrecSum(p,K)
%
%On return, either res is a faithful rounding of s:=sum(p), or
%
%   abs(res-s) < eps^K sum(abs(p)) .
%
%The default is K=2.
%
%The statements are also true in the presence of underflow. 
%The input vector p may be single or double precision.
%
%Implements Algorithm 4.1 from
%  S.M. Rump, T. Ogita, S. Oishi: Fast High Precision Summation,
%    Nonlinear Theory and Its Applications (NOLTA), IEICE, Vol. 1, No. 1, 2010.
%Requires (4L+3)n flops for L as computed in the algorithm. For double precsion
%and the default K=2, we have L=2 or L=3. Note that an implementation in C should
%use the modification as in Algorithm 5.2 in the cited paper to avoid extra memory.
%
%Reference implementation! Slow due to interpretation!
%

% written  09/29/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
% modified 01/24/12     S.M. Rump  comment
%

  if nargin==1
    K = 2;
  end
  
  if ~isreal(p)
    res = complex(PrecSum(real(p),K),PrecSum(imag(p),K));
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
  n = length(p);
  if n>nmax
    error(['maximum length of input vector ' int2str(nmax)])
  end
  
  if isa(p,'double'), prec='double'; else prec='single'; end
  u = 0.5*eps(prec);
  mu = sum(abs(p))/(1-2*n*u);
  if mu==0
    res = 0; 
    if rndold, setround(rndold); end
    return
  end
  Ms = NextPowerTwo(n+2);               % n+2 <= 2^M = Ms
  M = log2(Ms);
  L = ceil( ( K*log2(u) - 2 ) / ( log2(u) + M ) ) - 1;
  sigma = zeros(L+1,1);                 % sigma_k in the paper is sigma(k+1)
  tau = zeros(L,1);                     % hosts sum of leading terms
  sigma(1) = NextPowerTwo(mu);          % first extraction unit
  if ~isfinite(sigma)
    error('overflow occurred in PrecSum')
  end
  
  phi = Ms*u;                           % factor to decrease sigma
  for k=1:L                             % compute sigma(i)
    if sigma(k)>realmin(prec)
      sigma(k+1) = phi*sigma(k);
    else
      L = k-1; break
    end
  end
  
  if L==0                               % sigma_0 in underflow range
    res = sum(p);
    if rndold, setround(rndold); end
    return
  end
  
  for k=1:L                             % main loop: vertical version
    q = ( sigma(k) + p ) - sigma(k);    % [tau,p] = ExtractVector(sigma_{k-1},p);
    tau(k) = sum(q);                    % sum of leading terms
    p = p - q;                          % remaining terms
  end
  
  pi = tau(1); e = 0;
  for k=2:L                             % compute final result
    x = pi;
    pi = x + tau(k);                    % [pi,q]=FastTwoSum(pi,tau(k))
    q = tau(k) - ( pi - x );
    e = e + q;                          % fl-pt sum of errors
  end
  res = pi + ( e + sum(p) );            % final result
  
  if rndold
    setround(rndold)
  end
  