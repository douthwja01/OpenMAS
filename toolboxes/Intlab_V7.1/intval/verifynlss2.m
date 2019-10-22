function [ X , xs, K ] = verifynlss2(funfcn,xs,K,see,varargin)
%VERIFYNLSS2  Verified solution of nonlinear system for double and multiple roots
%
%The function comprises of various methods depending on the input parameters. Here
%  f:R->R of f:R^n->R^n is always a function to be called by f(x), possibly with
%  extra parameters, f(x,P1,P2,...). Optional input "see" is always to see intermediate
%  results.
%
%The univariate case f:R->R, xs is always an approximation to a K-fold root of f(x)=0.
%Optional output xs is an improved approximation.
%
%   [ X , xs , K ] = verifynlss2(f,xs,K,see,P1,P2,...)
%
%K>1:  There is an exactly K-fold root xhat of the perturbed equation
%        f(x) - e2*x^(K-2) - e3*x^(K-3) - ... - e(K-1)*x - e(K)
%      with xhat in X(1) and e(i) in X(i) for 2<=i<=K.
%
%K<-1: There are exactly |K| roots of the original equation f(x)=0 in the complex disc X.
%
%K=0:  The number of roots will be guessed, and then tried to include |K| roots of the original
%        equation f(x)=0 in the complex disc X. 
%
%Note that the disc is complex because a real double root of f may disappear into two
%complex roots by an arbitrarily small perturbation of f.
%For |K|=1 a unique simple root is included by verifynlss.
%
%
%The multivariate case f:R^n->R^n. Call:
%
%   [ XE , xs , K ] = verifynlss2(f,xs,K,see,P1,P2,...)
%
%Input is a function f:R^n->R^n, xs is an approximation xs to a double root of f(x)=0.
%
%optional input    K        equation to be regularized, 1<=K<=n
%                           []   default K=1
%                           -1   optimal K automatically determined
%                   see     see intermediate results 
%                   P1,...  extra parameters for function evaluation f(x,P1,P2,...)
%
%output             XE      inclusion of solution
%                   [xs;e]  improved approximation (column vector)
%                   K       choosen equation which was normalized
%
%There is a root xhat of the perturbed equation g(x)=0 with 
%   g_i(x) := f_i(x) for i~=K   and
%   g_K(x) := f_K(x)-ehat 
%such that g(xhat)=0 and Jacobian(g,xhat)=Jacobian(f,xhat) is singular.
%The inclusion XE satisfies xhat in XE(1:n) and ehat in XE(n+1).
%
%This means the function g nearby f has truely a double root. If the given 
%  function f has a double root, then ehat=0. However, the inclusion of ehat is
%  usually a narrow interval containing zero. Note that proving a function f to have
%  a double root is an ill-posed problem and outside the scope of verification methods.
%Starting with xs first some Newton steps are executed; if the approximation xs to
%  a double root is poor, the inclusion E may be not near zero.
%The result XE is NaN if no inclusion could be computed.
%
%Simple, one-dimensional nonlinear functions which can be written in one
%formula string, can be entered directly. For example,
%    X = verifynlss2(@(a)(sin(a)-1),1.5)
%evaluates an inclusion of the (double) zero of sin(a)-1=0 near as=1.5.
%
%See "help verifynlss" for references to the nonlinear system solver for simple roots.
%Inclusions for double roots of multivariate functions and multiple roots of
%univariate functions is based on
%  S.M. Rump, S. Graillat: Verified error bounds for multiple roots of systems 
%    of nonlinear equations, Numerical Algorithms, 54(3):359–377, 2009. 
%Inclusions of root cluster of univariate functions is based on
%  A. Neumaier: An Existence Test for Root Clusters and Multiple Roots, 
%    ZAMM, 68(6):256-257, 1988.
%and
%  S.M. Rump and S. Oishi: Verified Computation of a Disc containing exactly 
%    k roots of a univariate nonlinear function, Nonlinear Theory and Its
%    Applications (NOLTA), IEICE, 1(1), 1-8, 2009. (PDF, 108220 bytes)
%

%Multiple roots based on Taylor expansions:
% 
%                   k-1          k-1-nu
% g(x) := f(x) - sum     F_nu * X            with
%                   nu=1 
%         (k-1-j)              (k-j)                j-1  (k-1-nu)!         j-nu
% F_j := T       (xs) + (k-j)*T     (X)*(X-xs) - sum     --------- F_nu * X
%                                                   nu=1  (j-nu)!
%            (nu)
%It follows g    (xhat) = 0  for  0 <= nu <= k-1 .
%
                                                    
% written  04/14/09     S.M. Rump
% modified 08/02/12     S.M. Rump  comment
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  % store warning mode
  wng = warning;
  warning off
  
  % store standard function exception mode
  RealStdFctsExcptnMode = intvalinit('RealStdFctsExcptn',0);
  intvalinit('RealStdFctsExcptnNaN',0);
  
  if nargin>2
    if isempty(K)
      K = 2;
    end
  else
    K = 2;
  end

  xs = xs(:);
  if ( nargin<4 ) | isempty(see)
    see = 0;
  end
  
  % Convert to inline function as needed
  strfun = fcnchk(funfcn,length(varargin));

  if K==0                         % guess multiplicity
    k = zeros(1,3);
    for j=1:3
      y = feval(strfun,gradientinit(xs),varargin{:});
      xsnew = xs - (1:20)*y.x/y.dx;
      if ~( any(isinf(xsnew)) | any(isnan(xsnew)) );
        ynew = abs(feval(strfun,xsnew,varargin{:}));
        ynew(ynew==0)=realmin;
        [dummy,k(j)] = max(ynew(2:end)./ynew(1:end-1));
        if see
          disp(['guessed multiplicity ' int2str(k(j))'])
        end
        xs = xsnew(k(j));
      else
        if j==1
          error('guessing multiplicity failed')
        else
          k(j) = k(j-1);
        end
      end
    end
    if any(diff(k))
      warning('guessing multiplicity vague')
      if k(1)==k(2)
        K = -k(1);
      else
        K = -k(2);
      end
    else
      K = -k(1);
    end
  end
  
  % dimension
  n = length(xs);
  sqrt_n = sqrt(n);
  if ( n~=1 ) & ( K~=-1 ) & ( ( K<1 ) | ( K>n ) )
    error('invalid value of K: must be empty, -1, or between 1 and dimension')
  end
  
  if n==1                         % one nonlinear function

    % inclusion of root of (|K|-1)-st derivative
    [XX,xs] = verifynlss(funfcn,xs,abs(K)-1,see,varargin{:});
    if see
      disp(['inclusion of root of ' int2str(abs(K)-1) '-st derivative'])
      XX
    end
 
    if K>1                       % K-fold root of perturbed equation
    
      X = intval(zeros(K,1));
      X(1) = XX;
      xs = mid(XX);
      H = XX-xs;
      if K==2
        y = feval(strfun,intval(xs),varargin{:});
        Y = feval(strfun,gradientinit(XX),varargin{:});
        X(2) = y + Y.x*H;
      else
        y = feval(strfun,taylorinit(intval(xs),K-2),varargin{:});
        if see
          disp(['Taylor expansion at approximate solution xs=' num2str(xs)])
        end
        Y = feval(strfun,taylorinit(XX,K-1),varargin{:});
        X(2) = y{K-2} + (K-1)*Y{K-1}*H;
        for j=2:K-1
          nu = (1:j-1)';
          factor = intval(ones(j-1,1));
          for mu=nu
            factor(mu) = intval(prod(1:K-1-mu))/prod(1:j-mu);
          end
          X(j+1) = y{K-1-j} + (K-j)*Y{K-j}*H - sum(factor.*X(nu+1).*XX.^(j-nu));
        end
      end

    elseif K<-1              % |K| roots of original equation in complex disc X
      
      % Better to use f^(k)(xs)+f^(k+1)(X)(X-xs) rather than f^(k)(X) 
      % compute start interval
      K = -K;
      Y = cintval(XX);
      t = feval(strfun,taylorinit(intval(xs),K),varargin{:});     % Taylor xs
      T = feval(strfun,taylorinit(Y,K+1),varargin{:});            % Taylor Y
      gy = t{K}+abs((Y-XX)*T{K+1});
      if in(0,gy)
        X = intval(NaN);            % leading coefficient contains zero
      else
        P = polynom([gy 0 t{K-2:-1:0}]);                          % f^(k-1)(xhat)=0
        Z = XX + cintval(0,rootbound(P));                         % start interval
      
        % interval iteration
        k = 0; succ = 0;
        while ~succ & ( k<10 ) & ~( any(isnan(Z)) | any(isinf(Z)) )
          k = k+1;
          Y = Z*midrad(1,1e-15)+cintval(0,1e-324);                % epsilon-inflation
          t = feval(strfun,taylorinit(intval(xs),K),varargin{:}); % Taylor xs
          T = feval(strfun,taylorinit(Y,K+1),varargin{:});        % Taylor Y
          gy = t{K}+abs((Y-XX)*T{K+1});
          if in(0,gy)
            X = intval(NaN);        % leading coefficient contains zero
            break                   % stop iteration
          else
            P = polynom([gy 0 t{K-2:-1:0}]);                        % f^(k-1)(xhat)=0
            Z = XX + cintval(0,rootbound(P));                       % potential inclusion
            succ = in0(Z,Y);
            if see
              disp(['potential inclusion interval iteration'])
              disp(midrad(Z))
            end
          end
        end
        if succ
          X = Z;
          K = -K;
        else
          X = intval(NaN);
        end
      end
      
    else                            % |K|=1  verifynlss should have been used
      
      X = XX;
      
    end
    
  else                              % system of equations
    
    % Normalization component for kernel of Jacobian
    y = feval(strfun,gradientinit(xs),varargin{:});
    if K==-1                        % choose equation to be regularized
      try
        [L,U,p] = lu(mid(y.dx),'vector');
        [dummy,p] = sort(p);          % inverse permutation
      catch
        [L,U,p] = lu(mid(y.dx));
        [dummy,p] = max(p);
      end
      en = zeros(n,1); en(n) = 1;
      x0 = L(p,:)'\en;              % kernel vector:  A'*x0 ~ 0
      [dummy,K] = max(abs(x0));     % K-th component of kernel vector maximal
    end
    try
      [L,U,p] = lu(mid(y.dx)','vector');
      [dummy,p] = sort(p);            % inverse permutation
    catch
      [L,U,p] = lu(mid(y.dx)');
      [dummy,p] = max(p);
    end
    en = zeros(n,1); en(n) = 1;
    x0 = L(p,:)'\en;                % kernel vector:  A*x0 ~ 0
    [dummy,LL] = max(abs(x0));      % LL-th component of kernel vector maximal
    x0 = x0/x0(LL);                 % normalize LL-th component to 1

    % floating point Newton iteration for double root
    k = 0;
    xs = [ xs ; 0 ; x0(1:LL-1) ; x0(LL+1:n) ];  % add shift parameter and kernel vector
    xsold = xs;
    dxs = zeros(size(xs));
    dxsnew = abs(xs);
    me1 = zeros(n,1); me1(K) = -1;
    while ( any(dxsnew<.5*dxs) & ( norm(dxsnew)>sqrt_n*1e-14*norm(xs) ) & k<100 ) | ( k<3 )
      k = k+1;                    % at most 100, at least 3 iterations performed
      dxs = dxsnew;
      xsold = xs;
      y = mid(feval(strfun,hessianinit(xs(1:n)),varargin{:}));
      xs0 = [ xs(n+2:n+LL) ; 1 ; xs(n+LL+1:2*n) ];
      if K==1
        fxs = [ y.x(1)-xs(n+1) ; y.x(2:n) ; y.dx*xs0 ];
      else
        fxs = [ y.x(1:K-1) ; y.x(K)-xs(n+1) ; y.x(K+1:n) ; y.dx*xs0 ];
      end
      J = zeros(n);
      for i=1:n
        J(i,:) = xs0'*y(i).hx;
      end
      J = [ y.dx me1 zeros(n,n-1) ; J zeros(n,1) y.dx(:,1:LL-1) y.dx(:,LL+1:n) ];
      xs = xs - J\fxs;
      if see
        disp(['residual norm(f(xs_k)), floating point iteration ' sprintf('%d',k)])
        norm(y.x)
      end
      dxsnew = abs(xs-xsold);
    end
    
    % interval iteration
    R = inv(J);
    YY = feval(strfun,hessianinit(intval(xs(1:n))),varargin{:});
    xs0 = [xs(n+2:n+LL);1;xs(n+LL+1:2*n)];
    if K==1
      fxs = [ YY.x(1)-xs(n+1) ; YY.x(2:n) ; YY.dx*xs0 ];
    else
      fxs = [ YY.x(1:K-1) ; YY.x(K)-xs(n+1) ; YY.x(K+1:n) ; YY.dx*xs0 ];
    end
    Z = - R * fxs;
    X = Z;
    ready = 0; k = 0; kmax = 10;
    while ( ~ready ) & ( k<kmax ) & ( ~any(isnan(X)) )
      E = 0.1*rad(X)*hull(-1,1) + midrad(0,realmin);
      k = k+1;
      if see
        disp(['interval iteration ' sprintf('%d',k)])
      end
      Y = hull( X + E , 0 );       % epsilon inflation
      Yold = Y;
      YY = feval(strfun,hessianinit(xs(1:n)+Y(1:n)),varargin{:});
      J = intval(zeros(n));
      for i=1:n
        J(i,:) = (xs0+[Y(1:LL-1);0;Y(LL+1:n)])'*YY(i).hx;
      end
      J = [ YY.dx me1 zeros(n,n-1) ; J zeros(n,1) YY.dx(:,1:LL-1) YY.dx(:,LL+1:n) ];
      C = eye(2*n) - R * J;   % automatic gradients
      i=0;
      while ( ~ready ) & ( i<2 )   % improved interval iteration
        i = i+1;
        X = Z + C * Y;
        if any(isnan(X(:)))
          ready = 0;
          break
        end
        ready = all(all(in0(X,Y)));
        Y = intersect(X,Yold);
      end
    end
      
    if ready & isempty(find(isnan(Y) | isinf(Y) ))
      X = xs+Y;                    % verified inclusion
      X = X(1:n+1);
    else
      X = intval(repmat(NaN,n+1,1)); % inclusion failed
    end

    if nargin>1
      xs = xs(1:n+1);
    end

  end
    
  % restore warning and exception mode
  warning(wng)
  % restore out-of-range exception mode
  intvalinit(RealStdFctsExcptnMode,0);
  
  if rndold
    setround(rndold)
  end
  