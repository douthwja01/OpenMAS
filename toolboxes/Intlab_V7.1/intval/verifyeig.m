function [L,X] = VerifyEig(A,lambda,xs,B)
%VERIFYEIG      Verification of eigenvalue (cluster) near (lamda,xs)
%
%   [L,X] = VerifyEig(A,lambda,xs,B)
%
%Input: an approximate eigenvalue/eigenvector pair (lambda,xs).
%If B specified, generalized eigenproblem  A*x = lambda*B*x  is treated.
%On output, L contains exactly one eigenvalue and X contains one
%  eigenvector, which correspond to each other.
%
%For an eigenvalue cluster and/or multiple eigenvalue near lambda, 
%  xs(:,i), i=1:k, is an approximation to the corresponding invariant
%  subspace.
%On output, L contains (at least) k eigenvalues of A, and X
%  includes a base for the corresonding invariant subspace.
%By principle, clusters must be included in a complex interval L.
%
%Methods for simple eigenvalues are based on the solution of systems of nonlinear
%equations, see verifynlss, and have been presented in my Ph.D. thesis
%  S.M. Rump: Kleine Schranken für Matrixprobleme, Karlsruhe, 1980. 
%Methods for multiple eigenvalues are based on
%  S.M. Rump: Computational error bounds for multiple or nearly multiple eigenvalues, 
%    LAA 324, 209-226, 2001.
%

% written  07/15/99     S.M. Rump
% modified 10/20/01     S.M. Rump  Generalized eigenproblem added
% modified 06/26/02     S.M. Rump  output always of type intval, also for NaN
% modified 02/20/03     S.M. Rump  supress warnings
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
% modified 02/07/08     S.M. Rump  typo (thanks to Jiri Rohn)
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 09/04/11     S.M. Rump  comment in header
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  wng = warning;
  warning('off');

  [n k] = size(xs);

  [N,I] = sort(sum(abs(xs),2));
  u = I(1:n-k);
  v = I(n-k+1:n);
  midA = mid(A);

  if nargin==3                         % not generalized eigenproblem

    % one floating point iteration
    R = midA - lambda*speye(n);
    R(:,v) = -xs;
    y = R\(midA*xs-lambda*xs);
    xs(u,:) = xs(u,:) - y(u,:);
    lambda = lambda - sum(diag(y(v,:)))/k;

    R = midA - lambda*speye(n);
    R(:,v) = -xs;
    R = inv( R );
    C = A - intval(lambda)*speye(n);
    Z = - R * ( C * xs );
    C(:,v) = -xs;
    C = speye(n) - R * C;
    Y = Z;
    Eps = 0.1*mag(Y)*hull(-1,1) + midrad(0,realmin);
    m = 0;
    mmax = min( 15 * ( sum(sum(mag(Z(v,:))>.1)) + 1 ) , 20 );
    ready = 0;
    while ( ~ready ) & ( m<mmax ) & ( ~any(isnan(Y(:))) )
      m = m+1; 
      X = Y + Eps;
      XX = X;
      XX(v,:) = 0;
      Y = Z + C*X + R*(XX*X(v,:));
      ready = all(all(in0(Y,X)));
    end

  else                                 % generalized eigenproblem

    midB = mid(B);

    % one floating point iteration
    R = midA - lambda*midB;
    R(:,v) = -midB*xs;
    y = R\(midA*xs-lambda*midB*xs);
    xs(u,:) = xs(u,:) - y(u,:);
    lambda = lambda - sum(diag(y(v,:)))/k;

    R = midA - lambda*midB;
    R(:,v) = -midB*xs;
    R = inv( R );
    C = A - intval(lambda)*B;
    Z = - R * ( C * xs );
    C(:,v) = -B*xs;
    C = speye(n) - R * C;
    Y = Z;
    Eps = 0.1*mag(Y)*hull(-1,1) + midrad(0,realmin);
    m = 0;
    mmax = 15 * ( sum(sum(mag(Z(v,:))>.1)) + 1 );
    ready = 0;
    while ( ~ready ) & ( m<mmax ) & ( ~any(isnan(Y(:))) )
      m = m+1;
      X = Y + Eps;
      XX = X;
      XX(v,:) = 0;
      Y = Z + C*X + R*((B*XX*X(v,:)));
      ready = all(all(in0(Y,X)));
    end

  end

  if ready
    M = mag(Y(v,:));                         % eigenvalue correction
    if length(v)==1                          % one eigenvalue
      L = midrad(lambda,M);
    else                                     % eigenvalue cluster
      [Evec,Eval] = eig(M);
      [rho,index] = max(abs(diag(Eval)));
      Perronx = abs(Evec(:,index));
      setround(1);
      rad = max( ( M*Perronx ) ./ Perronx );   % upper bound for Perron root
      setround(0)
      L = cintval(midrad(lambda,rad));
    end
    Y(v,:) = 0;
    X = xs + Y;
  else
    % disp('no inclusion achieved')
    X = intval(NaN*ones(size(xs)));
    L = intval(NaN);
  end
  warning(wng)                         % reset warning mode
    
  if rndold
    setround(rndold)
  end
