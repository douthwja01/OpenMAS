function [res,err] = Dot_(varargin)
%DOT_         Dot products in K-fold precision
%
%   res = Dot_(A1,x1,...,An,xn,K)
%
%   [res,err] = Dot_(A1,x1,...,An,xn,K)
%
% res or res+err is sum(Ai*xi)
%     K optional, default K=2
%     K   >0  approximate result as if computed in K-fold precision
%         <0  inclusion of result as if computed in |K|-fold precision
%
% First factor Ai may be either matrix or scalar; in the latter case
%   res is of the size of xi.
% Second factor xi must be vector or scalar.
% All factors may be complex, but must be non-interval. 
% Sizes of every product Ai*xi must be the same.
%
% Precision K must be at least 2.
%
% Interval input makes only sense, if all intervals are degenerated (point intervals); 
%   therefore omitted.
%
% For speed and to fight interpretation overhead, 2 extra matrices are necessary; for
%   the product of to nxn matrices this would be n^3 memory, therefore this is not allowed.
%
% Examples are
%
%   res = Dot_(x',y);                  Approximate scalar product x'*y in quadruple precision
%   res = Dot_(A,x,-1,b);              Approximate residual A*x-b in quadruple precision
%   res = Dot_(A,x,-1,b,-2);           Inclusion of residual A*x-b in quadruple precision
%   res = Dot_(A,x1,A,x2,-1,b,3);      Approximate residual A*x1+A*x2-b in 3-fold precision
%
%For randomly generated full matrices of dimension n with b=a*randn(n,1) the following
%  table lists the computing times
%    t1   for ordinary residual  b - A*x
%    t2   with improved residual by lssresidual
%    t3   with quadruple precision residual by Dot_
%as well as the ratios to t1. All times are in seconds on a 800 MHz Pentium III Laptop.
%
%    n     t1      t2     t3    t2/t1  t3/t1     approximation of b-A*x
%--------------------------------------------
%   50   0.0002  0.001  0.003    3.1   14.6
%  100   0.0002  0.003  0.010   13.7   45.7
%  200   0.0005  0.010  0.045   20.5   89.5
%  500   0.0064  0.089  0.416   13.7   64.6
% 1000   0.0246  0.347  1.698   14.1   68.9
%
%The same table for verified inclusion of the residual b-A*x is as follows.
%
%    n     t1      t2     t3    t2/t1  t3/t1     inclusion of b-A*x
%--------------------------------------------
%   50   0.0012  0.001  0.004    0.9    3.0
%  100   0.0013  0.002  0.010    1.9    7.5
%  200   0.0020  0.012  0.048    6.1   23.8
%  500   0.0138  0.103  0.381    7.5   27.5
% 1000   0.0505  0.421  1.863    8.3   36.9
%
%Finally we display the achieved accuracy of the three methods for 
%  matrices of dimension 100 and different condition numbers, a randomly
%  generated right hand side b and computed solution x=A\b. 
%  We treat 100 test cases each and display the median and maximum componentwise
%  relative error of the results. Matrices are normed to 1.
%
%                A*x-b          lssresidual         Dot_
%condition  median maximum    median maximum    median maximum
%---------------------------------------------------------------
%    1e2    6.8e-16 1.0e+0   5.0e-21 8.9e-16   0.0e+0  4.7e-30 
%    1e5    1.8e-01 1.0e+0   1.8e-18 3.0e-01   5.0e-29 1.7e-27 
%   1e10    2.9e-01 1.0e+0   2.4e-06 6.8e-01   3.3e-24 1.8e-22 
%   1e14    2.8e-01 1.0e+0   2.1e-06 1.0e+0    2.7e-20 2.6e-18  
%
%The more we pay, the more we get.
%
% Implements algorithms Dot2 and DotK from
%   T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005
%

% written  11/08/03     S.M. Rump
% modified 11/30/03     S.M. Rump  complex data allowed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 08/08/04     S.M. Rump  sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/12/05     S.M. Rump  apapted to paper, improved performance
% modified 05/30/07     S.M. Rump  typo
% modified 08/07/10     S.M. Rump  upper case Dot_
% modified 06/06/11     S.M. Rump  result in two parts
% modified 08/26/12     S.M. Rump  global variables removed, rounding
% modified 09/05/12     S.M. Rump  sparse complex (thanks to M. Lange)
% modified 10/17/12     S.M. Rump  comments
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  % detect parameters
  len = length(varargin);
  if even(len)                      % default precision
    K = 2;
  else
    K = varargin{len};
  end
  idim = 0;
  
  % check precision K
  if abs(K)<2
    error('invalid parameter K for Dot_')
  end
  
  % detect dimension of result
  if prod(size(varargin{1}))==1
    [m n] = size(varargin{2});
  else
    m = size(varargin{1},1);        % m or n must be 1
    n = size(varargin{1},2);
  end
  
  % Initialization
  factor = 134217729;               % splitting factor 2^27+1
  if abs(K)==2
    res = zeros(m,1);
    if K>0
      s = 0;
    else
      sinf = 0;
      ssup = 0;
    end
  else
    res = [];
    S = [];
  end
  
  % Compute array to sum up
  for i=1:floor(nargin/2)
    
    % compute next product A*x
    A = varargin{2*i-1};
    x = varargin{2*i};
    sA = size(A);
    sx = size(x);
    if sx(2)~=1
      error('second factor in Dot_ must be vector')
    end
    
    if ( prod(sA)==1 ) & ( sx(1)~=m )
      error('dimensions of summands in Dot_ do not match')
    end

    if ( prod(sA)==1 ) & ( abs(A)==1 )   % first factor A scalar +/- 1
      
      % TwoProduct(A,x)
      if A==1
        [res,q] = sumK([res x],2);
      else
        [res,q] = sumK([res -x],2);
      end
      if abs(K)==2
        if K>0
          s = s + q;
        else
          setround(-1)
          sinf = sinf + q;
          setround(1)
          ssup = ssup + q;
          setround(0)
        end
      else
        S = [S q];
      end
      idim = idim + 1;
      
    else                                    % first factor non-scalar or ~=  +/- 1
      
      % matrix-vector product or dot product
      if ( prod(sA)~=1 ) & ( sA(2)~=sx(1) )
        error('inner dimensions in Dot_ do not match')
      end
      if issparse(A)
        [iA,jA,ssA] = find(A);
      end
      C = factor*A;
      A1 = C - ( C - A );                       % upper part of A
      A2 = A - A1;                              % A = A1+A2 exact splitting
      
      C = factor*x;
      x1 = C - ( C - x );                       % upper part of x
      x2 = x - x1;                              % x = x1+x2 exact splitting
      
      I = sqrt(-1);
      if sA(1)~=1                               % m~=1: matrix-vector product
        if issparse(A)                          % factors s.t. A*x = sum( (A1+A2).*(x1+x2) , 2 )
          x1 = sparse(iA,jA,x1(jA),m,sA(2));
          x2 = sparse(iA,jA,x2(jA),m,sA(2));
        else
          x1 = x1(:,ones(1,m)).';
          x2 = x2(:,ones(1,m)).';
        end
        if isreal(A) | isreal(x)                % A or x real
          if issparse(A)
            h = sparse(iA,jA,ssA.*x(jA),m,sA(2));
          else
            h = A .* ( x(:,ones(1,m)).' );
          end
          idim = idim + sx(1);
        else                                    % both A and x complex
          if issparse(A)
            h = [ sparse(iA,jA,ssA.*real(x(jA)),m,sA(2)) , sparse(iA,jA,I*(ssA.*imag(x(jA))),m,sA(2)) ];
          else
            xx = x(:,ones(1,m)).';
            h = [ A.*real(xx) , A.*(I*imag(xx)) ];
          end
          idim = idim + 2*sx(1);
        end
      elseif sA(2)~=1
        if isreal(A) | isreal(x)                % A or x real
          h = A .* x.';
          x1 = x1.';
          x2 = x2.';
          idim = idim + sx(1);
        else                                    % both A and x complex
          h = [ A.*real(x).' , A.*(I*imag(x)).' ];
          x1 = x1.';
          x2 = x2.';
          idim = idim + 2*sx(1);
        end
      else
        if isreal(A) | isreal(x)                % A or x real
          h = A * x;
          idim = idim + 1;
        else                                    % both A and x complex
          h = [ A*real(x) , A*(I*imag(x)) ];
          idim = idim + 2;
        end
      end
      
      % TwoProduct(A,x): A*x = h + r
      if isreal(A) | isreal(x)                % A or x real
        r = A2.*x2 - ((( h - A1.*x1 ) - A2.*x1 ) - A1.*x2 );
      else                                    % both A and x complex
        r = [A2.*real(x2),A2.*(I*imag(x2))] - ((( h - [A1.*real(x1),A1.*(I*imag(x1))] ) - ...
          [A2.*real(x1),A2.*(I*imag(x1))] ) - [A1.*real(x2),A1.*(I*imag(x2))] );
      end
      
      if issparse(h)                    % compress h
        h = compress(h);
      end
      if issparse(r)                    % compress r
        r = compress(r);
      end
      
      % Summation vector and partial summation for |K|=2
      if abs(K)==2
        [res,q] = sumK([res h],2);      % error-free
        sq2 = size(q,2);
        sr2 = size(r,2);
        if sq2>sr2
          r = [r zeros(size(r,1),sq2-sr2)];
        elseif sq2>sr2
          q = [q zeros(size(q,1),sr2-sq2)];
        end
        if K>0
          s = s + sum(q+r,2);
        else
          setround(-1)
          sinf = sinf + sum(q+r,2);
          setround(1)
          ssup = ssup + sum(q+r,2);
          setround(0)
        end
      else
        [res,q] = sumK([res h],2);
        S = [S q r];
      end
      
    end
    
  end

  % Compute result in K-fold precision
  if abs(K)==2                          % quadruple precision
    if K>0
      res = res + s;
    else
      INTLAB_INTVAL_ETA = realmin*eps;  % smallest positive denormalized fl-pt
      setround(-1)
      resinf = sinf - 5*idim*INTLAB_INTVAL_ETA;
      setround(1)
      ressup = ssup + 5*idim*INTLAB_INTVAL_ETA;
      err = infsup(resinf,ressup);
      if nargout==1
        res = res + err;
      end
    end
  else                                  % |K|-fold precision
    if K>0
      res = Sum_([S res],K,2);
    else
      INTLAB_INTVAL_ETA = realmin*eps;  % smallest positive denormalized fl-pt
      [res,q] = sumK([S res],-K);       % error-free transformation
      setround(-1)
      resinf = sum(q,2) - 5*idim*INTLAB_INTVAL_ETA;
      setround(1)
      ressup = sum(q,2) + 5*idim*INTLAB_INTVAL_ETA;
      err = infsup(resinf,ressup);
      if nargout==1
        res = res + err;
      end
    end
  end
    
  if rndold
    setround(rndold)
  end

  
function [res,q] = sumK(a,K)
%SUM          Summation of "a" in K-fold precision along 2nd dimension
%             res is approximate sum, sum of errors in q
%             res + sum(q,2) = sum(a,2) in exact arithmetic

  % number of summands
  N = size(a,2);
  
  % summation (i.e. quadruple for K==2)
  for i=1:K-1
    pi = cumsum(a,2);
    z = diff(pi,1,2);
    q = ( pi(:,1:N-1) - ( pi(:,2:N) - z ) ) + ( a(:,2:N) - z );
    res = pi(:,N);
    if i~=K-1
      a = [q res];
    end
  end

function a = compress(a)
%COMPRESS     horizontal compress of array "a" such that sum(a,2) remains unchanged

  [i,j,s]=find(a.'); 
  rows = sum(spones(a),2);
  m = max(rows);
  if m==0
    a = sparse([],[],[],size(a,1),m);
    return
  end
  rows = nonzeros(rows);
  index = cumsum([1;rows(1:end-1)]);
  h = ones(size(i));
  h(index(2:end))=-rows(1:end-1)+1;
  a = sparse(j,cumsum(h),s,size(a,1),m);