function [L,X] = structeig(As,lambda,xs)
%STRUCTEIG      Verification of eigencluster near (lamda,xs) for structured input matrix
%
%   [L,X] = StructEig(A,lambda,xs)
%
%For an eigenvalue cluster near lambda, where xs(:,i), i=1:k is an
%  approximation to the corresponding invariant subspace
%
%Matrices may be treated subject to different structures producing different 
%  results. Try, for example, inclusions of the double eigenvalue zero in the
%  following example:
%
%    e = 1e-3;
%    A = midrad( toeplitz([0 1+e -e/2 1+e]),1e-4); 
%    [v,d] = eig(A.mid); xs = v(:,2:3); lambda = d(2,2);
%    X1 = verifyeig(A,lambda,xs);
%    X2 = structeig(structure(A,'symmetric'),lambda,xs);
%    X3 = structeig(structure(A,'symmetricToeplitz'),lambda,xs);
%    X4 = structeig(structure(A,'generalToeplitz'),lambda,xs);
%    X5 = structeig(structure(A,'persymmetric'),lambda,xs);
%    X6 = structeig(structure(A,'circulant'),lambda,xs);
%    res = [ X1 X2 X3 X4 X5 X6 ]', rad(res)
%
%On output, L contains (at least) k eigenvalues of A, and X
%  includes a base for the corresonding invariant subspace
%By principle, L is a complex interval.
%
%Methods for simple eigenvalues are based on the solution of systems of nonlinear
%equations, see verifynlss, and have been presented in my Ph.D. thesis
%  S.M. Rump: Kleine Schranken fE Matrixprobleme, Karlsruhe, 1980,
%multiple eigen values based on
%  S.M. Rump: Computational error bounds for multiple or nearly
%    multiple eigenvalues, LAA, 324:209-226, 2001,
%and structured input data based on
%  S.M. Rump: Rigorous Sensitivity Analysis for Systems of Linear and
%    Nonlinear Equations. Math. Comput., 54(10):721-736, 1990,
%see also
%  S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic.
%    Acta Numerica, 19:287-449, 2010. 
%

% written  07/15/99     S.M. Rump
%

  if ~isstruct(As)
    error('input matrix must be structured')
  end

  [n k] = size(xs);
  rndold = getround;                   % get current rounding mode
  setround(0)

  [N,I] = sort(sum(abs(xs),2));
  u = I(1:n-k);
  v = I(n-k+1:n);
  A = reshape( mid(As.Phi*As.p) ,n,n);
  midA = mid(A);

  % one floating point iteration
  R = midA - lambda*speye(n);
  R(:,v) = -xs;
  y = R\(midA*xs-lambda*xs);
  xs(u,:) = xs(u,:) - y(u,:);
  lambda = lambda - sum(diag(y(v,:)))/k;

  R = midA - lambda*speye(n);
  R(:,v) = -xs;
  R = inv( R );
  [i,j,s] = find(As.Phi);
  q = floor((i-1)/n);
  kk = size(As.Phi,2);
  Z = ( R*intval(xs) ) * lambda;
  for nu=1:k
    Z(:,nu) = Z(:,nu) - ...
      ( intval(R) * sparse(i-n*q,j,s.*xs(q+1,nu),n,kk) ) * As.p ;
  end
  C = A - intval(lambda)*speye(n);

  C(:,v) = -xs;
  C = speye(n) - R * C;
  Y = Z;
  Eps = 0.1*abs(Y)*hull(-1,1) + midrad(0,realmin);
  m = 0;
  mmax = 15 * ( sum(sum(abs(Z(v,:))>.1)) + 1 );
  ready = 0;
  while ( ~ready ) & ( m<mmax ) & ( ~any(isnan(Y(:))) )
    m = m+1;
    X = Y + Eps;
    XX = X;
    XX(v,:) = 0;
    Y = Z + C*X + R*(XX*X(v,:));
    ready = all(all(in0(Y,X)));
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
    disp('no inclusion achieved')
    X = NaN*ones(size(xs));
    L = NaN;
  end
  setround(rndold)                     % reset old rounding mode
  