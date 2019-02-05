function res = norm(a,r)
%NORM         Inclusion of norm of a matrix or vector
%
%   res = norm(a,r)
%
%   r in {1,2,inf}         for vector a
%   r in {1,2,inf,'fro'}   for matrix a
%
% default for r is 2
%
%Implementation for the spectral norm of real or complex matrices based on
%
%   S.M. Rump: Verified bounds for singular values, in particular for the
%      spectral norm of a matrix and its inverse, BIT, 51(2):367-384, 2011.
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged, redesign
% modified 01/18/06     S.M. Rump  Matlab bug
% modified 02/17/10     S.M. Rump  matrix spectral norm
% modified 03/01/10     S.M. Rump  isspd
% modified 03/16/10     S.M. Rump  symmetry/Hermitian ensured
% modified 11/25/10     S.M. Rump  intersection [0,inf] (thanks to 
%                                    Helmut Podhaisky, Halle)
% modified 03/24/11     S.M. Rump  approximation improved (thanks to Y. Watanabe, 
%                                    Fukuoka for pointing to this problem)
% modified 01/20/12     S.M. Rump  abs(A) in normposA 
%
  if nargin==1
    r = 2;
  end
  [m,n] = size(a);
  
if isequal(r,'inf'), r=inf; end
  if ( m==1 ) | ( n==1 )        % vector a
    if r==1
      res = sum(abs(a));
    elseif r==inf
      res = abs(a);
      res = intval(max(inf(res)),max(sup(res)),'infsup');
    elseif ( r==2 ) | isequal(r,'fro')
      %VVVV a = a(:);           % column vector
      s.type = '()'; s.subs = {':'}; a = subsref(a,s);
      %AAAA Matlab V5.2 bug fix
      res = sqrt(abs(sum(a'*a)));  % use abs to make sure sum is nonnegative
    else                        % invalid parameter for vector norm
      error('invalid request for norm')
    end
  else                          % matrix a
    if r==1                     % column sum norm
      suma = sum(abs(a),1);
      res = intval(max(inf(suma)),max(sup(suma)),'infsup');
    elseif r==inf               % row sum norm
      suma = sum(abs(a),2);
      res = intval(max(inf(suma)),max(sup(suma)),'infsup');
    elseif isequal(r,'fro')     % Frobenius norm
      %VVVV a = a(:);           % staggered columns
      s.type = '()'; s.subs = {':'}; a = subsref(a,s);
      %AAAA Matlab V5.2 bug fix
      res = sqrt(abs(sum(a'*a)));  % use abs to make sure sum is nonnegative
      
    elseif r==2                 % matrix spectral norm
      nfull = 50;               % full eigendesomposition up to this dimension
      acc = 1e-6;               % output accuracy

      e = 1e-30;                % in comments: input matrix A
      if 1+e==1-e               % fast check for rounding to nearest
        rndold = 0;
      else
        rndold = getround;
        setround(0)
      end
      
      % split interval matrix into midpoint and radius
      if a.complex              % complex input
        if a.rad==0             % point matrix
          normR = intval(0);
        else
          normR = normposA(a.rad,acc);  % norm(A.rad)<=normR.sup
        end
        a = a.mid;              % careful: fl-pt matrix A.mid
      else                      % real input
        if isequal(a.inf,a.sup) % point matrix
          normR = intval(0);
          a = a.inf;            % careful: fl-pt matrix A.mid
        else
          setround(1)
          amid = a.inf + 0.5*(a.sup-a.inf);
          a = amid-a.inf;       % radius matrix
          setround(0)
          normR = normposA(a,acc);  % norm(A.rad)<=normR.sup
          a = amid;             % careful: fl-pt matrix A.mid
        end
      end
      
      % treat midpoint matrix a (may be sparse)      
      if isequal(a,a')          % symmetric or Hermitian, use Corollary 3.2
        symm = 1;
      else                      % general matrix, use Theorem 3.3
        symm = 0;
      end
      
      if n<=nfull               % compute full decomposition
        a = full(a);
        if symm                 % symmetric/Hermitian A.mid
          [V,D] = eig(a);       % approximate decomposition of A.mid
          alphaV = norm(eye(n)-V'*intval(V),1);
          [dummy,k] = max(abs(diag(D)));    % spectral radius in D_kk
          D = V'*(a*intval(V)); % "almost" similiarity transformation of A.mid
          N = max(mag(gershgorin(D)))/(1-alphaV);  % upper bound for norm(D)
          x = intval(V(:,k));
          N = infsup(inf(norm(a*x)/norm(x)),N.sup); % bounds for norm(A.mid)
          res = N + intval(-normR.sup,normR.sup,'infsup');  % bounds for norm(A)
        else                    % general A.mid
          [V,D] = eig(a'*a);    % approximate decomposition of A.mid'*A.mid
          alphaV = norm(eye(n)-V'*intval(V),1);
          [dummy,k] = max(abs(diag(D)));    % spectral radius of A.mid'*A.mid in D_kk
          x = intval(V(:,k));
          D = a*intval(V);
          D = D'*D;             % "almost" similarity transformation of A.mid'*A.mid
          N = sqrt(max(mag(gershgorin(D)))/(1-alphaV));   % upper bound for norm(A.mid)
          N = infsup(inf(sqrt(norm(a'*(a*x))/norm(x))),N.sup);  % bounds for norm(A.mid)
          res = N + intval(-normR.sup,normR.sup,'infsup');      % bounds for norm(A) 
        end
        if res.inf<0            % intersection with [0,inf]
          res.inf = 0;
        end
        if rndold
          setround(rndold)
        end
        return
      end
      
      % larger dimension n>nfull
      if ~symm                  % prepare for general case
        a = a'*intval(a);       % A.mid'*A.mid
        a = tril(a)+tril(a,-1)';  % make sure result is symmetric/Hermitian
        if a.complex            % complex matrix
          normR = normR + normposA(a.rad,acc);
          a = a.mid;            % midpoint of A.mid'*A.mid
        else                    % real matrix
          setround(1)
          amid = a.inf + 0.5*(a.sup-a.inf);
          a = amid-a.inf;
          setround(0)
          normR = normR + normposA(a,acc);
          a = amid;             % midpoint of A.mid'*A.mid
        end 
      end
      
      if issparse(a)
        p = symamd(a);          % reordering to reduce fill-in
        a = a(p,p);
      end

      % fl-pt iteration to approximate norm(A.mid) or norm(A.mid'*A.mid)
      x = sum(abs(a),1)';
      normest = norm(x);
      for i=0:1                     % two trial vectors to avoid numerical artefacts
        if i, x = x.*(1+.2*randn(n,1)); end
        normold = 0;
        while abs(normest-normold)>acc*normest/20 % improved accuracy
            normold = normest;
            x = x/normest;
            x = a*x;
            normest = norm(x);
        end
        if ~i, normest0 = normest; end
      end
      normest = max(normest,normest0);
      
      % verification of upper bound
      phi = 20;                 % relax accuracy
      alpha = (1+phi*acc)*normest;   % alpha >~ norm(A.mid) or norm(A.mid)^2
      setround(-1)
      b = alpha*speye(n)-a;
      setround(0)
      if symm                   % verify alpha*I+/-A.mid s.p.d.
        while 1
          if isspd(b,0,0)
            setround(-1)
            b = alpha*speye(n)+a;
            setround(0);
            if isspd(b,0,0), break, end
          end
          phi = 2*phi;          % relax bound
          alpha = (1+phi*acc)*alpha; % iteration: increase alpha (never happens)
          setround(-1)
          b = alpha*speye(n)-a;
          setround(0)
        end
      else                      % verify alpha*I-A.mid'*A.mid s.p.d.
        while ~isspd(b,0,0)
          phi = 2*phi;          % relax bound
          alpha = (1+phi*acc)*alpha; % iteration: increase alpha (never happens)
          setround(-1)
          b = alpha*speye(n)-a;
          setround(0)
        end
      end
      
      % verified lower bound
      x = intval(x);
      alpha = hull( norm(a*x)/norm(x) , alpha );  % lower bound on norm(a)
      if ~symm              % make sure norm(a) <= alpha for original a
        alpha = sqrt(alpha);    % lower bound on norm(A.mid)
      end
      normR = intval(-normR.sup,normR.sup,'infsup');  % make sure normR is symmetric
      res = alpha + normR;
      if res.inf<0              % intersection with [0,inf]
        res.inf = 0;
      end
      if rndold                 % reset rounding
        setround(rndold)
      end
      
    else                        % invalid parameter for matrix norm
      error('invalid request for norm')
    end
  end
  
  
function normbnd = normposA(A,normrelerr)
% upper bound for norm of nonnegative matrix A up to normrelerr according to (3.1/2)
  normbnd = inf;
  cont = 1;
  x = sum(A,1)';
  if ~any(x)                % take care of A==0
    normbnd = 0;
    return
  end
  setround(1)
  while cont & ( normbnd~=0 )   % use Collatz' lemma
    normold = normbnd;
    y = A'*(A*x);           % y_i=0 only for x==0
    normbnd = max(y./x);    % zero column of A produces NaN in y./x
    x = y/normbnd;
    cont = ( abs(normbnd-normold)>normrelerr*normbnd );
  end
  setround(0)
  normbnd = sqrt(intval(normbnd));
