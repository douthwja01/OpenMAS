function [ X , xs , k ] = verifynlss(funfcn,xs,param,see,varargin)
%VERIFYNLSS   Verified solution of nonlinear system
%
%   [ X , xs , k ] = verifynlss(f,xs,param,see,P1,P2,...)
%
%f is name of function, to be called by f(xs), xs is approximation.
%If param is 'g' or 's', then f:R^n->R^n ;
%if param is 'h', then f:R^n->R.
%Result X is an inclusion of xhat near xs with f(xhat)=0 or f'(xhat)=0.
%  Result is NaN if no inclusion could be computed.
%By default, the function is expanded by gradients. If param is 
%  specified 's', slopes are used instead.
%If param is specified 'h', then hessians are used and f'(xhat)=0
%  is solved for f:R^n->R.
%This function is designed for simple roots; for multiple roots use
%  "verifynlss2".
%
% optional input    param   'g'  use gradient, proves uniqueness
%                           's'  use slopes, better, but w/o uniqueness
%                           'h'  solve f'(x)=0
%                           'gSparseSPD', 'sSparseSPD', 'hSparseSPD'
%                               same as above but using sparse linear
%                               system solver
%                           'gSparse', 'sSparse', 'hSparse' same as above but
%                               using sparse linear system solver to normal
%                               equations, thus limiting cond(Jacobian) to about 1e8
%                           k  integer, solve k-th derivative f^(k)(x)=0 (only for
%                                univariate functions)
%                   see     see intermediate results 
%                   P1,...  extra parameters for function evaluation
%                           f(x,P1,P2,...)
%
% optional output   xs      improved approximation (column vector)
%                   k       interval iteration steps
%
%
%Simple, one-dimensional nonlinear functions which can be written in one
%formula string, can be entered directly. The unknown must be 'x'. E.g.,
%    X = verifynlss('x*exp(x)-1',.6)
%evaluates an inclusion of the zero of x*exp(x)-1=0 near xs=.6.
%
%Nonlinear system solver based on the Krawcyzk operator, see
%  R. Krawczyk: Newton-Algorithmen zur Bestimmung von Nullstellen mit
%    Fehlerschranken, Computing 4, 187-201, 1969.
%  R.E. Moore: A Test for Existence of Solutions for Non-Linear Systems,
%    SIAM J. Numer. Anal. 4, 611-615, 1977.
%with modifications for enclosing the error with respect to an approximate
%solution, an iteration scheme, epsilon-inflation and an improved interval
%iteration as in
%  S.M. Rump: Solving Algebraic Systems with High Accuracy, in "A New
%    Approach to Scientific Computation", eds. U. Kulisch and W. Miranker,
%    Academic Press, 51-120, 1983.
%  S.M. Rump: Verification methods for dense and sparse systems of equations, 
%    in : J. Herzberger (ed.), Topics in Validated Computations - Studies in 
%    Computational Mathematics, Elsevier, Amsterdam, 63-136, 1994.
%
%Using gradient verifies existence and uniqueness of a zero of the nonlinear
%  function within the inclusion interval. This also implies multiplicity 1 of
%  the zero, and nonsingularity of the Jacobian at the zero.
%Using slopes implies existence but not uniqueness of a zero. This allows
%  inclusion of zero clusters and multiple zeros. For details, see
%    S.M. Rump: Inclusion of Zeros of Nowhere Differentiable n-dimensional
%      Functions, Reliable Computing, 3:5-16 (1997).
%Using param='h' uses hessians, i.e. gradients for f'.
%

% written  10/16/98     S.M. Rump
% modified 10/12/99     S.M. Rump  output NaN-vector in case of failure,
%                                  interval iteration stops if NaN occurred
% modified 06/26/02     S.M. Rump  output always of type intval, also for NaN
% modified 01/09/04     S.M. Rump  stopping criterion for fl-pt iteration changed following
%                                    a proposal by Arnold Neumaier
% modified 04/04/04     S.M. Rump  hessians and zeros of f' added, improved fl-pt stopping criterion
%                                    set round to nearest for safety
%                                    sparse Jacobian/Hessian
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/16/07     S.M. Rump  see intermediate result for parameter 'h'
% modified 09/06/07     S.M. Rump  warning supressed
% modified 10/27/07     S.M. Rump  check for NaN in iteration
% modified 11/14/07     S.M. Rump  solve for normal equations added
% modified 03/09/08     S.M. Rump  check for parameter corrected
% modified 09/23/08     S.M. Rump  stdfct exception corrected
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/19/08     S.M. Rump  StdFctsException ignore/NaN, check for inf
% modified 05/24/09     S.M. Rump  k-th derivative using Taylor package
% modified 10/05/09     S.M. Rump  parameter check
% modified 08/02/12     S.M. Rump  comment verifynlss2
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

  xs = xs(:);
  if ( nargin<3 ) | isempty(param)
    param = 'g';
  else
    if ~ischar(param)
      if isnumeric(param) & (length(param)==1) & ( param>=0 ) & ( param==round(param) )
        kderiv = param;
        if kderiv==0
          param = 'g';
        elseif kderiv==1
          param = 'h';
        else
          if length(xs)>1
            error('roots of higher derivatives only for univariate functions')
          end
          param = 'k';
        end
      else
        error('third parameter of verifynlss must be char or integer')
      end
    else
      param = lower(param);
    end
  end
  
  if ( nargin<4 ) | isempty(see)
    see = 0;
  end
  
  % Convert to inline function as needed
  strfun = fcnchk(funfcn,length(varargin));
  
  % floating point Newton iteration
  n = length(xs);
  sqrt_n = sqrt(n);
  xsold = xs;
  dxs = zeros(size(xs));
  dxsnew = abs(xs);
  k = 0;
  while ( any(dxsnew<.5*dxs) & ( norm(dxsnew)>=sqrt_n*1e-14*norm(xs) ) & k<100 ) | ( k<3 )
    k = k+1;                     % at most 100, at least 3 iterations performed
    dxs = dxsnew;
    xsold = xs;
    if isequal(param(1),'g') | isequal(param(1),'s') 
      y = feval(strfun,gradientinit(xs),varargin{:});
      xs = xs - y.dx\y.x;
    elseif isequal(param(1),'h')
      y = feval(strfun,hessianinit(xs),varargin{:});
      xs = xs - y.hx\y.dx';
    else
      y = feval(strfun,taylorinit(xs,kderiv+1),varargin{:});
      xs = xs - (kderiv+1)*y{kderiv+1}\y{kderiv};
    end
    if see
      if isequal(param(1),'g') | isequal(param(1),'s')     % solution of f(x)=0
        disp(['residual norm(f(xs_k)), floating point iteration ' sprintf('%d',k)])
        norm(y.x)
      elseif isequal(param(1),'h') % solution of f'(x)=0
        disp(['residual norm(f''(xs_k)), floating point iteration ' sprintf('%d',k)])
        norm(y.dx)
      else                         % solution of f^(kderiv)(x)=0
        disp(['residual norm(f^(' int2str(kderiv) ')(xs_k)), floating point iteration ' sprintf('%d',k)])
        norm(y{kderiv})
      end
    end
    dxsnew = abs(xs-xsold);
  end
  
  % check full or sparse case
  if length(param)>1
    sparse_ = ~isempty(findstr(param,'sparse'));
  else
    sparse_ = 0;
  end
  
  if sparse_
    
    % sparse interval iteration
    if isequal(param(1),'h')
      y = feval(strfun,hessianinit(intval(xs)),varargin{:});
      Z = - y.dx';
    else
      Z = - feval(strfun,intval(xs),varargin{:});
    end
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
      if isequal(param(1),'g')        % use gradients
        x = gradientinit(xs+Y);
        y = feval(strfun,x,varargin{:});
        M = y.dx;                  % automatic gradients
      elseif isequal(param(1),'s')    % use slopes
        x = slopeinit(xs,xs+Y);
        y = feval(strfun,x,varargin{:});
        M = y.s;                   % automatic slopes
      else                         % use hessians
        x = hessianinit(xs+Y);
        y = feval(strfun,x,varargin{:});
        M = y.hx;                  % automatic hessians
      end
      i=0;
      solve_spd = isequal(param(end-2:end),'spd');
      while ( ~ready ) & ( i<2 )   % improved interval iteration
        i = i+1;
        if solve_spd
          % solve s.p.d. system
          X = verifylss(M.mid,Z+mag(M.rad)*abs(Y));
        else
          % solve normal equations
          X = verifylss(intval(M.mid')*M.mid,M.mid'*(Z+mag(M.rad)*abs(Y)));
        end
        ready = all(all(in0(X,Y)));
        Y = intersect(X,Yold);     % intersection, therefore M still valid
      end
    end
    
  else
    
    % full interval iteration
    if isequal(param(1),'g')  | isequal(param(1),'s')   % f(x)=0
      R = inv(y.dx);
      Z = - R * feval(strfun,intval(xs),varargin{:});
    elseif isequal(param(1),'h')      % f'(x)=0
      R = inv(y.hx);
      y = feval(strfun,hessianinit(intval(xs)),varargin{:});
      Z = - R * y.dx';
    else                              % f^(kderiv)(x)/k!=0
      R = 1/((kderiv+1)*y{kderiv+1});
      Y = feval(strfun,taylorinit(intval(xs),kderiv),varargin{:});
      Z = - R * Y{kderiv};
    end
    X = Z;
    ready = 0; k = 0; kmax = 10;
    while ( ~ready ) & ( k<kmax ) & ( ~any(isnan(X)) )
      E = 0.1*rad(X)*hull(-1,1) + midrad(0,realmin);
      k = k+1;
      if see
        disp(['interval iteration ' sprintf('%d',k)])
      end
      Y = hull( X + E , 0 );          % epsilon inflation
      Yold = Y;
      if isequal(param(1),'g')
        x = gradientinit(xs+Y);
        y = feval(strfun,x,varargin{:});
        C = eye(n) - R * y.dx;        % automatic gradients
      elseif isequal(param(1),'s')    % use slopes
        x = slopeinit(xs,xs+Y);
        y = feval(strfun,x,varargin{:});
        C = eye(n) - R * y.s;         % automatic slopes
      elseif isequal(param(1),'h')    % use hessian
        x = hessianinit(xs+Y);
        y = feval(strfun,x,varargin{:});
        C = eye(n) - R * y.hx;        % automatic hessians
      else                            % f^(kderiv)(x)/k!=0              
        x = taylorinit(xs+Y,kderiv+1);
        y = feval(strfun,x,varargin{:});
        C = 1 - R * ( (kderiv+1)*y{kderiv+1} );   % automatic Taylor expansion
      end
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
    
  end
  
  if ready & isempty(find(isnan(Y) | isinf(Y) ))
    X = xs+Y;                    % verified inclusion
  else
    X = intval(repmat(NaN,n,1)); % inclusion failed
  end
  
  % restore warning and exception mode
  warning(wng)
  % restore out-of-range exception mode
  intvalinit(RealStdFctsExcptnMode,0);
  
  if rndold
    setround(rndold)
  end
  