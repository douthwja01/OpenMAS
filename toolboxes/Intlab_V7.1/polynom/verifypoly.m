function [ X , k , exact ] = verifypoly(P,xs,k)
%VERIFYPOLY   Verified inclusion of polynomial root
%
%General algorithm to include roots of (real or complex polynomials). Call
%
%   [ X , k , exact ] = verifypoly(P,xs,k)
%
%Input   P      polynomial
%        xs     approximation to root(s) of P
%        k      (optional) number of roots near xs
%Output  X      interval containing at least k roots of P
%        k      if input k specified: 
%                   usually the same value except extraordinary situations (with warning)
%               if input k not specified:   determined by the routine
%        exact  true if exactly k roots of P in X
%
%The algorithm should always deliver an including interval X except leading coefficient 
%  containing zero.
%
%Methods based on
%  S.M. Rump: Ten methods to bound multiple roots of polynomials,  
%    J. Comput. Appl. Math. (JCAM), 156:403Å-432, 2003.
%and the literature cited over there, especially
%  P. Batra: Abschatzungen und Iterationsverfahren fÅE Polynom-Nullstellen, Ph.D. Thesis,
%    Institute for Reliable Computing, Hamburg University of Technology, 1999.
%  A. Neumaier: Enclosing clusters of zeros of polynomials,
%    J. Comput. Appl. Math. (JCAM), 156:389-401, 2003.
%

% written  09/21/02     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  midP = mid(P);
  n = degree(midP);
  wngold = warning;
  warning off;
  
  if nargin==2                                        % guess multiplicity
    z = roots(midP);
    absP = polynom(abs(vector(midP)));
    absPabsxs = absP{abs(xs)};
    k = 1;
    for m=1:n
      P_m = polynom(binom(n:-1:m,m).*midP(n:-1:m));   % m-th derivative of P divided by m!
      r = ( eps * absPabsxs/abs(P_m{xs}) ) ^ (1/m);
      ms = sum( abs(z-xs)<2*r );
      if ms==m
        k = m;
        sens = r;
        break
      end
    end
  else
    if k~=1
      P_k = polynom(binom(n:-1:k,k).*midP(n:-1:k));   % m-th derivative of P divided by m!
      absP = polynom(abs(vector(midP)));
      sens = ( eps * absP{abs(xs)}/abs(P_k{xs}) ) ^ (1/k);
    end
  end    
  
  if k==1                                             % inclusion method for simple zero
    
    % floating point Newton iteration
    midPs = midP';
    if nargin==2                                      % initialize with nearest zero
      [dummy,index] = min(abs(z-xs));
      xs = z(index(1));
    end
    xsold = xs;
    i = 0;
    while ( norm(xs-xsold)>1e-10*norm(xs) & i<10 ) | i<1
      i = i+1;                                        % at most 10, at least 1 iteration performed
      xsold = xs;
      midP_xs = midP{xs};
      midPs_xs = midPs{xs};
      corr = midP_xs/midPs_xs;
      if abs(corr)>0.5*abs(midP_xs)
        corr = 0.5*sign(corr)*abs(midP_xs);
      end
      xs = xs - corr;
    end
    
    % interval iteration
    Ps = intval(P)';
    R = 1/midPs_xs;
    Z = - R * P{intval(xs)};
    X = Z;
    E = 0.1*rad(X)*hull(-1,1) + midrad(0,realmin);
    ready = 0; j = 0; jmax = 10;
    while ( ~ready ) & ( j<jmax ) & ( ~any(isnan(X)) )
      j = j+1;
      Y = hull( X + E , 0 );                          % epsilon inflation
      Yold = Y;                      
      C = 1 - R * Ps{xs+Y};
      i=0;    
      while ( ~ready ) & ( i<2 )                      % improved interval iteration
        i = i+1;
        X = Z + C * Y;
        ready = in0(X,Y);
        Y = intersect(X,Yold);
      end
    end
    if ready
      exact = 1;
      X = xs+Y;                                       % verified inclusion
      [X,kk] = finalcheck(n,P,X);
      warning(wngold)
      if ( nargin==3 ) & ( kk~=0 ) & ( k~=kk )
        s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
        warning(s)
      end
      if kk~=0, k = kk; end  
      setround(rndold)
      return
    else                                              % no bound for simple zero, calculate a priori bound
      exact = 0;
      Xs = intval(xs);
      absPXs = intval(mag(P{Xs}));
      r = sup( absPXs^(1/intval(n)) );
      r2 = sup(n*absPXs/abs(Ps{Xs}));
      if ~isnan(r2)
        r = min(r,r2);
      end
      X = midrad(xs,r);                               % at least 1 zero in X
      [X,kk] = finalcheck(n,P,X);
      warning(wngold)
      if ( nargin==3 ) & ( kk~=0 ) & ( k~=kk )
        s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
        warning(s)
      end
      if kk~=0, k = kk; end  
      setround(rndold)
      return
    end
    
  end
  
  
  % inclusion method for multiple zero
  
  if nargin==3                                        % zero approximations not yet computed
    z = roots(midP);
  end
  [dummy,index] = sort(abs( z-xs ));
  z = z(index);
  c = mean(z(1:k));
  Q = pshift(P,intval(c));
  
  if in(0,Q(k))                                       % extreme numerical situation or wrong k, ue backup
    [X,kk,exact] = backup(n,k,Q,P,z,c,sens);          % at least or exactly kk zeros in X
    warning(wngold)
    if ( nargin==3 ) & ( k~=kk )
      s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
      warning(s)
    end
    k = kk;  
    setround(rndold)
    return
  end
  
  if all(Q(0:k-1)==0)                                 % exact k-fold root
    exact = 1;
    X = intval(c);
    warning(wngold)  
    setround(rndold)
    return
  end
  
  % compute approximation to radius (lower bound) by Cauchy polynomial
  q = polynom([mig(Q(k)) -mag(Q(k-1:-1:0))]);
  v = abs(q(k-1:-1:0)/q(k));
  v(k) = abs(q(0)/(2*q(k)));
  r = 2 * max( v.^(1./(1:k)) );   % Fujiwara root bound for q, off by maximal factor 2
  rold = 2*r;
  qs = q';
  acc = 1e-2;
  while ( r>0 ) & ( r<100*sens) & ( abs(r-rold)>acc*abs(r+rold) ) 
    rold = r;
    qr = q{r};
    r = r - qr/qs{r};                                 % approximates smallest positive root of Q
  end
  r1 = r;
  
  % one (damped) Newton correction on r to produce upper bound for radius
  qq = polynom([mag(Q(n:-1:k+1)) -mig(Q(k)) mag(Q(k-1:-1:0))]); 
  qqs = qq';
  qqr1 = qq{r};
  r = r - min(max(-r,qqr1/qqs{r}),0);                 % correction at least zero, at most r
  
  % check radius
  if r<0
    [X,kk,exact] = backup(n,k,Q,P,z,c,sens);          % no better bound possible, at least or exactly k zeros in X
    warning(wngold)
    if ( nargin==3 ) & ( k~=kk )
      s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
      warning(s)
    end
    k = kk;  
    setround(rndold)
    return
  end
  
  % Data for Pellet's test
  QQ = polynom([abs(Q(n:-1:k+1)) -abs(Q(k)) abs(Q(k-1:-1:0))]); 
  % sup(QQ{r})<0  implies  at exactly k zeros in D(c,r)
  
  % check approximation r
  success = ( sup(QQ{r})<0 );                         % Pellet's test
  
  if ~success                                         % initial guess not valid, the rare case
    
    % approximate root of qq'
    qqss = qqs';
    rold = 2*r;
    acc = 1e-2;
    i = 0;
    while ( abs(r-rold)>acc*abs(r+rold) ) & ( i<5 ) & ( r<100*sens )
      i = i+1;
      rold = r; 
      r = r - qqs{r}/qqss{r};
    end
    
    % one secant step
    qqr = qq{r};
    r = r - qqr/(qqr-qqr1)*(r-r1);
    
    % check radius
    if r<0
      [X,kk,exact] = backup(n,k,Q,P,z,c,sens);        % no better bound possible, at least or exactly k zeros in X
      warning(wngold)
      if ( nargin==3 ) & ( k~=kk )
        s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
        warning(s)
      end
      k = kk;  
      setround(rndold)
      return
    end
    
    % check approximation r
    success = ( sup(QQ{r})<0 );                       % Pellet's test
    
    if ~success                                       % no bounds by this procedure, calculate a priori bound
      [X,kk,exact] = backup(n,k,Q,P,z,c,sens);        % at least or exactly k zeros in X
      warning(wngold)
      if ( nargin==3 ) & ( k~=kk )
        s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
        warning(s)
      end
      k = kk;  
      setround(rndold)
      return
    end
    
  end
  
  % improvement of bound
  exact = 1;
  rho = 2*r;
  acc = 1e-5;
  imax = 5; i = 0;
  while ( abs(r-rho)>acc*abs(r+rho) ) & success & ( i<imax )
    i = i+1;
    rho = r;
    r = r - qq{r}/qqs{r};                             % approximates smallest positive root of q
    if r<0
      X = midrad(c,rho);
      warning(wngold)  
      setround(rndold)
      return
    end
    success = ( sup(QQ{r})<0 );                       % Pellet's test
  end
  
  X = midrad(c,rho);
  if rho>10*sens                                     % bound possibly improvable, use backups
    [Y,kk,exactY] = backup(n,k,Q,P,z,c,sens);         % at least or exactly k zeros in Y
    warning(wngold)
    if rad(Y)<rho
      if ( nargin==3 ) & ( k~=kk )
        s = ['verifypoly: ' int2str(kk) ' roots near approximation detected, but ' int2str(k) ' specified.'];
        warning(s)
      end
      k = kk;
      X = Y;
      exact = exactY;  
      setround(rndold)
      return
    end
  end
  warning(wngold)
    
  setround(rndold)

  
  function [X,k,exact] = backup(n,k,Q,P,z,c,sens)
  % backup function, compute upper bound for circle contaning at least k roots of P
  
  if ~in(0,Q(k))
    bound = 1;
    % compute van Vleck bound
    QvV = polynom([mig(Q(k)) mag(Q(k-1:-1:0))]);   % van Vleck polynomial
    setround(1)
    for i=0:k-1
      QvV(i) = -( QvV(i)*sup(binom(intval(n-i),k-i)) );
    end
    v = abs(QvV(k-1:-1:0)/QvV(k));
    v(k) = abs(QvV(0)/(2*QvV(k)));
    r = 2 * max( v.^(1./(1:k)) );   % Fujiwara root bound for QvV, off by at most factor 2
    rold = 2*r;
    QvVs = QvV';
    acc = 1e-2;
    imax = 5; i = 0;
    while ( abs(r-rold)>acc*abs(r+rold) ) & ( i<imax )
      i = i+1;
      rold = r;
      QvV_r = QvV{r};
      setround(-1)
      QvVs_r = QvVs{r};
      delta = QvV_r/QvVs_r;
      setround(1)
      r = r - delta;                                  % upper bound for root
    end
    X = midrad(c,r);
    exact = 0;
  else
    bound = 0;
  end
  
  if ( ~bound ) | ( rad(X)>2*sens ) 
    % van Vleck bound failed or too coarse, compute bound by Neumaier's modification of Gershgorin circles
    
    % compute vector p used for partial fraction expansion
    P_n = P(n);
    q = ones(n,1);
    Z = intval(z);
    for i=1:n
      qq = Z-z(i);
      qq(i) = 1;
      q = q.*qq;
    end
    p = P{Z}./q;
    
    % Neumaier's improved Gershgorin circles
    M = z - n/2*p/P_n;                            % midpoints (interval)
    R = abs(n/2*p/P_n);                           % Radii (interval)
    D = cintval(M + midrad(0,sup(R)));            % discs D_j according to Neumaier (8)
    MM = abs( M*ones(1,n) - ones(n,1)*M.' );      % midpoint distance matrix
    MR = R*ones(1,n) + ones(n,1)*R';              % radii sum matrix
    I = ones(n,1)*(1:n);
    I( inf(MM) > sup(MR) ) = 0;                   % I(i,j)~=0 => intersect(D_i,D_j) empty
    
    % compute cluster index set C of Neumaier-Gershgorin circles
    [dummy,index] = min(mag(z-intval(c)));
    C = setdiff(unique(I(index,:)),0);
    Cold = {};
    while ~isequal(C,Cold)
      Cold = C;
      C = setdiff(unique(I(C,:)),0);
    end
    lenC = length(C);
    
    mY = mean(mid(M(C))); 
    Y = cintval( midrad( mY , max(mag(M(C)-mY+R(C))) ) );
    lambda = inf;
    for j=C(:)'
      index = 1:n;
      index(C) = [];
      lambda = min( lambda , inf(sum(real( sum(p(index)./(D(j)-z(index)))/P_n ))) );
    end
    beta = 1+lambda;
    if beta*n/lenC>1
      d = p(C)/(2*beta)*lenC/P_n;
      MC = z(C)-d;
      midY = mean(mid(MC));
      YY = cintval( midrad( midY , max(mag(MC-midY+d)) ) );
      if YY.rad<Y.rad
        Y = YY;
      end
    end
    
    if bound
      if rad(Y)<rad(X)
        X = Y;
        exact = 1;
        k = length(C);
      end
    else
      X = Y;
      exact = 1;
      k = length(C);
    end
    
    success = ~isnan(X);
    c = mid(X);
    R = intval(rad(X));
    while success                                 % improve radius by Neumaier-Rouche, useful for large n
      R = R/2;
      d = c-intval(z);
      d2 = abs(d).^2;
      t = mig( real( 1 + sum( (conj(d).*p) ./ ( d2 - R^2 ) ) ) );
      s = mag( p ./ ( d2 - R^2 ) );
      setround(1)
      s = sup( R * sum(s) );
      setround(0)
      success = ( t>s );
      if success
        X = midrad(c,sup(R));
      end
    end
    
  end
  
  [X,kk] = finalcheck(n,P,X);
  if kk~=0
    k = n;
    exact = 1;
  end
    
  
  function [X,k] = finalcheck(n,P,X)
  % check the extremes: final bound should be better than root bound
  % k~=0  <=>  bound replaced
  
  k = 0;
  P_n = P(n);
  if isa(P_n,'intval')
    if in(0,P(n))                               % no bound possible
      return
    end
  else
    if abs(P_n)<1e-308                          % no bound
      return
    end
  end
  
  Pabs = polynom([mig(P(n)) -mag(P(n-1:-1:0))]);
  setround(1)
  v = abs(Pabs(n-1:-1:0)/Pabs(n));
  v(n) = abs(Pabs(0)/(2*Pabs(n)));
  r = 2 * max( mag(v.^(1./intval(1:n))) );     % Fujiwara root bound for Pabs
  if r<2*rad(X)                                 % apply some Newton iterations
    Pmag = Pabs';
    acc = 1e-2; 
    rold = 2*r;
    imax = 3; i = 0; 
    while ( abs(r-rold)>acc*abs(r+rold) ) & ( i<imax )
      i = i+1;
      rold = r;
      r = r - inf(Pabs{intval(r)}/Pmag{intval(r)});
    end
    setround(0)
    if r<rad(X)
      X = midrad(0,r);
      k = n;
    end
  else
    setround(0)
  end
