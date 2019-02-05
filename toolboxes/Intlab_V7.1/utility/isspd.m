function res = isspd(A,accurate,mmd)
%ISSPD        logical function: Matrix A is positive definite
%
%   res = isspd(A)
%
%Given a real symmetric or complex Hermitian matrix A, 
%
%   res   1  Matrix A is proved to positive definite
%         0  positive definiteness could not be verified
%
%If A is an interval matrix, then the above holds for every symmetric
%  (Hermitian) matrix contained in A
%
%The verification is performed based on a floating-point Cholesky
%decomposition, see 
%  S.M. Rump: Verification of positive definiteness, 
%    BIT Numerical Mathematics, 46:433-452, 2006. 
%
%This is fast, but for ill-conditioned matrices verification may not
%succeed. The call
%
%   res = isspd(A,1)
%
%is slower but may be successful for ill-conditioned matrices. It is
%based on
%
%  S.M. Rump: Validated Solution of Large Linear Systems, Computing
%    Supplementum 9, 191-212, 1993
%and
%  S.M. Rump: Verification of Positive Definiteness. BIT Numerical 
%    Mathematics, 46:433-452, 2006.
%

%Extra parameter:
%  mmd   0  input A already sorted w.r.t. minimum degree symamd
%        1  input A not presorted (default)
%

% written  01/21/05     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 01/18/06     S.M. Rump  dimension check
% modified 05/09/07     S.M. Rump  norm(Arad) improved, thanks to J. Rohn
%                                  symmetry check
% modified 10/06/08     S.M. Rump  diag(A).inf used, thanks to L. Kolev
%                                  diagonal check and dimension 1
% modified 09/08/12     S.M. Rump  literature
%

  % symmetry check
  if ~isequal(A,A')
    error('input matrix not symmetric/Hermitian')
  end

  % diagonal check
  d = diag(A); 
  if any( mid(d)<=0 )               % negative elements on diagonal
    res = 0;
    return
  end
  d = mig(d);                       % isspd(A_) => isspd(A) where A_ is A with
  if any( d<=0 ) | any( mid(d)<=0 ) % diagonal replaced by diag(A).inf 
    res = 0;
    return
  end
  
  % constants
  n = dim(A);
  if n==1
    res = 1;                        % A is positive real scalar, otherwise
    return                          % diagonal check would have failed
  end
  Eps = 2^(-53);
  Eta = 2^(-1074);
  if nargin==1
    accurate = 0;
  end
  if nargin<=2
    mmd = 0;
  end
  
  e = 1e-30;
  if 1+e==1-e                       % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  if isintval(A)                    % treat interval input
    Arad = rad(A);
    A = mid(A);
    A(1:n+1:n^2) = d;               % replace diag(A) by diag(A).inf with radius 0
    Arad(1:n+1:n^2) = 0;
  else
    M = 0;
  end
  
  % scaling
  maxdiagA = max(d);                % d must be positive
  mindiagA = min(d);
  if ( maxdiagA/mindiagA>sqrt(n) ) & ~( ( maxdiagA>1e100 ) | ( mindiagA<1e-100 ) )
    % apply van der Sluis scaling
    wng = warning;
    warning off
    ver = version;
    if isequal(ver(1:3),'5.3') & ( ~isreal(d) )  % Matlab 5.3 doesn't like log2 of
      d = 2.^(-ceil(0.5*log(d)/log(2)));         % diagonal of complex matrix
    else
      d = 2.^(-ceil(0.5*log2(d)));
    end
    warning(wng)
    D = spdiags( d ,0,n,n );          % D_ii are powers of 2
    A = D*A*D;                        % 0.25 <= abs(A_ii) < 1
    maxdiagA = 1;
    if exist('Arad','var')
      Arad = D*Arad*D;
    end
  end
  
  % Minimum degree sorting
  if mmd
    if exist('symamd','file')
      p = symamd(A);
    else
      if n<2000                       % symmmd can be very slow
        p = symmmd(A);
      end
    end
    if exist('p','var')               % sorting performed
      A = A(p,p);
      if exist('Arad','var')
        Arad = Arad(p,p);
      end
    end
  end
  
  if exist('Arad','var')
    % approximation of norm(rad(A)) = rho(rad(A))
    x = ones(n,1);
    m = 1;
    M = 2;
    iter = 0;
    setround(1)             % rounding upwards, Arad>=0, x>=0, use Collatz' bound
    while 1
      iter = iter+1;
      y = Arad*x;
      index = ( y==0 );
      if any(index)         % rows/cols can be removed
        Arad(index,:) = [];
        Arad(:,index) = [];
        x(index) = [];
        y(index) = [];
        m = 0;              % take care of zero components
      else
        m = inf;
      end
      if ~isempty(y)        % x and y entirely nonzero        
        x = y./x;
        m = min(m,min(x));  % take care of zero components
        M = max(x);
        x = y/max(y);
      else
        M = 0;
        break
      end
      if ( abs((M-m)/(m+M))<.1 ) | ( iter>10 )
        break
      end
    end
    % bound of norm(rad(A)) = rho(rad(A)) <= M
    setround(0)
  end
  
  if accurate                         % explicit computation of residual
    
    % approximation of smallest singular value and approximate solution
    % for positive definite A
    [C,p] = chol(A);
    if p==0                       % Cholesky decomposition of A.mid succeeded
      % approximation s of smallest eigenvalue
      y = C\(sum(abs(C),1)');     % initial approximation from normest
      s = n/(y'*y);
      k=0; notready=1;
      while notready | ( k<3 )    % approximation s for sigma_min(A.mid)
        k = k+1;
        sold = s;
        z = ( C'\y )*s;
        y = C\z;
        s = (z'*z)/(y'*y);
        notready = abs(s-sold)>.0001*abs(s+sold);
      end
      if s>M        % it may work
        s = .8*(s-M);
        setround(-1)
        A = A - s*speye(n);
        setround(0)
        [C,p] = chol(A);
        if p==0           % Cholesky decomposition of shifted matrix succeeded
          res = ( s>norm(C'*intval(C)-A) );
        else              % Cholesky decomposition of shifted matrix failed
          res = 0;
        end
      else
        res = 0;                  % s-M not positive   
      end
    else                          % Cholesky decomposition failed
      res = 0;
    end

  else                                % proof based on pure fl-pt compuation
    
    [i,j] = find(A);
    index = find(diff(j));
    t = [0 ; (2:n)'-i(index+1)];        % max #elts left of (or above) diag(A)
    if any( t > 1e15 )                  % dimension check
      res = 0;
      setround(rndold)
      return
    end

    % hull of A
    setround(1)                         % set rounding upwards
    if n<67108861                       % alpha_k/(1-alpha_k) < (k+1)*Eps
      alpha = (t+3)*Eps;
    else
      alpha = (t+2)*Eps;
      alpha = alpha./(1-alpha);
      alpha = alpha./(1-alpha);
    end
    d = sqrt(alpha.*diag(A));

    % Upper bound for norm(dA) and shift
    c = d'*d + ( M + 6*n*(n+maxdiagA)*Eta );

    % floating point Cholesky
    setround(-1)
    A = A - c*speye(n);

    setround(0)
    [R,p] = chol(A);
    res = ( p==0 );

  end
  
  setround(rndold)
