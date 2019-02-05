function resN = NearSum(p)
%NEARSUM      Computes the rounding to nearest result of sum(p)
%
%   resN = NearSum(p)
%
%On return, resN is sum(p) rounded to nearest, also in the presence
%  of underflow. Input vector p may be single or double precision.
%
%Implements Algorithm 7.3 from
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
    resN = complex(NearSum(real(p)),NearSum(imag(p)));
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
    error(['maximum length of input vector for NearSum ' int2str(nmax) '.'])
  end

  eta = succ(0*p(1));
  [tau1,tau2,p,sigma,Ms] = Transform(p,0);
  tau2s = tau2 + sum(p);
  [res,delta] = FastTwoSum(tau1,tau2s);     % res+delta = tau1+tau2s
  if delta==0                               % res = fl(s)
    resN = res;
    if rndold, setround(rndold); end
    return
  end
  R = tau2 - ( res - tau1 );                % s-res = R+sum(p)
  if delta<0                                % fl(s) in {pred(res),res}
    gamma = pred(res)-res;                  % res+gamma=pred(res)
    if gamma==-eta                          % s = res
      resN = res;
      if rndold, setround(rndold); end
      return
    end
    deltas = gamma/2;                       % mu := res+deltas = M^-(res)
    deltass = TransformK(p,R-deltas,sigma,Ms);  % s-mu = R-deltas+sum(p)
    if deltass>0                            % s > M^-(res)
      resN = res;
    elseif deltass<0                        % s < M^-(res)
      resN = pred(res);
    else                                    % s = M^-(res)
      resN = res + deltas;
    end
  else                                      % fl(s) in {res,succ(res)}
    gamma = succ(res)-res;                  % res+gamma=succ(res)
    if gamma==eta                           % s = res
      resN = res;
      if rndold, setround(rndold); end
      return
    end
    deltas = gamma/2;                       % mu := res+deltas = M^+(res)
    deltass = TransformK(p,R-deltas,sigma,Ms);  % s-mu = R-deltas+sum(p)
    if deltass>0                            % s > M^+(res)
      resN = succ(res);
    elseif deltass<0                        % s < M^+(res)
      resN = res;
    else                                    % s = M^+(res)
      resN = res + deltas;
    end
  end
  
  if rndold
    setround(rndold)
  end
  