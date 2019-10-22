function [I,k] = verifyquad(funfcn,a,b,tol,see,mmax,varargin)
%QUAD         Verified quadrature using Romberg's scheme - rudimentary implementation
%
%The input function must be real and at least 4 times differentiable. The quadrature
%  does NOT work if f or its derivatives have poles in or near the interval [a,b].
%  Not much optimization of the parameters is done; if the routine works well, there
%  is at least a verified inclusion of the integral available.
%
%   [ I , k ] = verifyquad(f,a,b,tol,see,mmax,varargin)
%
%The univariate real function f is integrated from a to b; must at least be 4 times differentiable.
%
% optional input    tol     anticipated (relative) tolerance (default 1e-6)
%                   see     see intermediate results 
%                   mmax    maximal order of Romberg scheme <=16; note that derivatives 
%                             up to 2*mmax+2 are computed (default mmax=4)
%                   P1,...  extra parameters for function evaluation
%                           f(x,P1,P2,...)
% optional output   k       number of subintervals
%
%A simple example:
%
%   Q = verifyquad('exp(x*sin(x))',0,10)
%
%Routine verifyquad is by no means optimized, just written straightforwardly. But sometimes
%it is even faster (and much more accurate) than the Matlab built-in function quad:
% f = @(x)(sinh(exp(x))), a = 0; b = 5; 
%   tic, Approx = quad(f,a,b), toc, 
%   tic, Incl = verifyquad(f,a,b), toc
%produces
% f = 
%     @(x)(sinh(exp(x)))
% Warning: Maximum function count exceeded; singularity likely.
% > In quad at 100
% Approx =
%   3.3124e+032
% Elapsed time is 0.281967 seconds.
% intval Incl = 
%   1.0e+061 *
%     9.6709
% Elapsed time is 0.161995 seconds.
%
%Note that verifyquad is almost twice as fast, and the approximate value is off
%  by several orders of magnitude. But a warning is given.
%It may also happen that no warning is given:
% f = @(x)(sin(x+exp(x))), a = 0; b = 8; 
%   tic, Approx = quad(f,a,b), toc, 
%   tic, Incl = verifyquad(f,a,b), toc
% f = 
%     @(x)(sin(x+exp(x)))
% Approx =
%     0.2511
% Elapsed time is 0.178671 seconds.
% intval Incl = 
%     0.3474
% Elapsed time is 1.125084 seconds.
%
%Note, however, that the approximate value is only correct to one figure, and no warning is given.
%
%Based on Romberg's rule with error term by
%   F. L. Bauer, H. Rutishauser, and E. Stiefel: New Aspects in Numerical Quadrature, 
%   Proc. Symp. Appl. Math. Vol. XV, AMS 1963, pp. 198-218.
%

% written  05/29/09     S.M. Rump
% modified 11/30/09     S.M. Rump  function check and vectorize
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  % store standard function exception mode
  RealStdFctsExcptnMode = intvalinit('RealStdFctsExcptn',0);
  intvalinit('RealStdFctsExcptnNaN',0);

  % Convert to inline function as needed
  try
    if isa(funfcn,'function_handle')
      strfun = fcnchk(eval(vectorize(funfcn)),length(varargin));
    else
      strfun = fcnchk(vectorize(funfcn),length(varargin));
    end
  catch
    strfun = fcnchk(funfcn,length(varargin));
  end

  if ( nargin<4 ) | isempty(tol)
    tol = 1e-6;                         % default (absolute) tolerance 1e-6
  end
  tol = max(tol,1e-16);                 % avoid tiny tolerances
  D = intval(b)-a;                      % inclusion of b-a
  if ( nargin<5 ) | isempty(see)
    see = 0;
  end
  if ( nargin<6 ) | isempty(mmax)
    mmax = 4;                           % up to 2*mmax+2 derivatives
  end
  if mmax>16
    warning('maximum order of Romberg scheme is mmax=16.')
    mmax = 16;
  end
 
  % store warning mode
  wng = warning;
  warning off
  
  I = infsup(-inf,inf);
  relerrold = inf;
  kmin = max(5,mmax);                   % start with at least 2^kmin subintervals
  
  % bounds for taylor^(2*m+2) for m=1..mmax
  xi = linspace(a,b,17);
  TF_ = feval(strfun,taylorinit(infsup(xi(1:end-1),xi(2:end)),2*mmax+2),varargin{:});
  if any(isnan(TF_.t(:))) | any(isinf(TF_.t(:)))  % second try with many subintervals
    xi = linspace(a,b,1025);
    TF_ = feval(strfun,taylorinit(infsup(xi(1:end-1),xi(2:end)),2*mmax+2),varargin{:});
  end
  if any(isnan(TF_.t(:))) | any(isinf(TF_.t(:)))  % derivative calculation failed
    X = feval(strfun,xi,varargin{:});
    if any(isnan(X(:))) | any(isinf(X(:)))
      error('function evaluation failed over interval [a,b]')
    else
      error('evaluation of higher derivatives failed, try smaller maximum order')
    end
  end
  for m=1:mmax
    TF{m} = infsup(min(TF_{2*m+2}.inf),max(TF_{2*m+2}.sup));
  end
  if see
    disp(['Inclusions of ' int2str(mmax) '-th derivative'])
    TF{mmax}
  end
  
  % first T{kmin..kmin-mmax+1}(m) for m=1
  k = kmin+1;
  h = 2^(-k)*D;                         % inclusion of (b-a)/2^k
  Xi = a + ( (0:2^k)*(2^(-k)) ) * D;    % inclusion of grid points a+i*(b-a)/2^k
  w = [ 1 repmat([4 2],1,2^(k-1)) ];    % no rounding errors
  w(end) = 1;
  y = feval(strfun,Xi,varargin{:})/3;
  T{k-1}(1) = h * sum(w.*y);
  for j=1:m-1
    w = w(1:(2^(k-j)+1));
    w(end) = 1;
    h = 2*h;                            % inclusion of (b-a)/2^k
    T{k-j-1}(1) = h * sum(w.*y(1:2^j:end));
  end
  
  % create tableaux up to mmax
  for m=2:mmax
    for k=kmin-mmax+1:kmin-m+1
      T{k}(m) = ( 4^m*T{k+1}(m-1) - T{k}(m-1) ) / ( 4^m-1 );
    end
  end
  
  % first error term for m=1    [ B(m) = +/- Bernoulli(2*m+2) ]
  B = intval([1 1 1 5 691 7 3617 43867 174611 854513 236364091 8553103 ...
                 23749461029 8615841276005 7709321041217 2577687858367]) ./ ...
             [30 42 30 66 2730 6 510 798 330 138 2730 6 870 14322 510 6];
  m = mmax;
  k = kmin-m+1;
  E = 2^(-2*k*(m+1)-m*(m+1)) * TF{m}*D^(2*m+3);
  err = B(m)*E;
  I = T{k}(m) - err;
  if see
    disp('First inclusion')
    I
  end
  
  % initialize Romberg iteration
  if in(0,I)
    relerrI = diam(I);
  else
    relerrI = relerr(I);
  end
  j = 0;
  
  while ( relerrI>tol ) & ( relerrI<relerrold ) & ( kmin+j<20 )
    j = j+1;
    Iold = I;
    relerrold = relerrI;
    k = kmin+j+1;
    % create entry T{kmin+j}(1)
    h = 2^(-k)*D;                         % inclusion of (b-a)/2^k
    Xi = a + ( (0:2^k)*(2^(-k)) ) * D;    % inclusion of grid points a+i*(b-a)/2^k
    w = [ 1 repmat([4 2],1,2^(k-1)) ];    % no rounding errors
    w(end) = 1;
    y = feval(strfun,Xi,varargin{:})/3;
    k = k-1;
    T{k}(1) = h * sum(w.*y);
    % update last row of Romberg table
    for m=2:mmax
      T{k-m+1}(m) = ( 4^m*T{k-m+2}(m-1) - T{k-m+1}(m-1) ) / ( 4^m-1 );
    end
    % new inclusion
    m = mmax;
    k = k-m+1;
    E = 2^(-2*k*(m+1)-m*(m+1)) * TF{m}*D^(2*m+3);
    err = B(m)*E;
    I = T{k}(m) - err;
    if see
      disp('New inclusion')
      I
    end
    if in(0,I)
      relerrI = diam(I);
    else
      relerrI = relerr(I);
    end
  end
  
  if relerrI>relerrold
    I = Iold;
    k = k-1;
  end
 
  % restore warning and exception mode
  warning(wng)
  % restore out-of-range exception mode
  intvalinit(RealStdFctsExcptnMode,0);
  
  if rndold
    setround(rndold)
  end
  