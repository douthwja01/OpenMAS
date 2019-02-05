function res = Sum_(a,K,dim)
%SUM_         Summation in K-fold precision
%
%   res = Sum_(A,K,dim)
%
% input  A    array to be summed, must be non-interval, up to 2 dimensions
%        K    >0   approximate result as if computed in K-fold precision 
%                    (default K=2, i.e. quadruple precision)
%             <0   inclusion of true sum, as if computed in |K|-fold precision
%        dim  sum along dimension dim (optional, default 1)
% output res  result of summation
%
% Precision K must be at least 2.
%
% Interval input makes only sense, if all intervals are degenerated (point intervals); 
%   therefore omitted.
%
% Functionality as Matlab function sum. For accuracy and timing see routine Dot_.
%
% Implements algorithms Sum2 and SumK from
%   T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%     SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005
%

% written  11/08/03     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/08/05     S.M. Rump  apapted to paper
% modified 08/07/10     S.M. Rump  upper case Dot_
%

  [m n] = size(a);
  
  % determine precision of computation
  if nargin==1
    K = 2;
  end
  
  % determine dimension to sum along
  if nargin<=2
    if m==1
      dim = 2;
    else
      dim = 1;
    end
  end
  
  % default K=2
  if isempty(K)
    K = 2;
  end
  
  % check precision K
  if abs(K)<2
    error('invalid parameter K for Sum_')
  end
  
  % treat complex input
  if ~isreal(a)
    res = Sum_(real(a),K,dim) + sqrt(-1)*Sum_(imag(a),K,dim);
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  % determine number of summands and possible transposition
  if dim==1
    N = m;
    transp = 0;
  else
    N = n;
    a = a.';
    transp = 1;
  end
  
  if ( m==1 ) | ( n==1 )    % vector sum
    if  ( ( ( m==1 ) & ( dim==1 ) ) | ( ( n==1 ) & ( dim==2 ) ) )
      if transp
        res = a.';
      else
        res = a;
      end
      return
    end
  end
  
  % summation along dimension 1
  for i=1:abs(K)-2          % K-2 error-free vector transformations
    pi = cumsum(a);
    z = diff(pi);
    a = [ ( pi(1:N-1,:) - ( pi(2:N,:) - z ) ) + ( a(2:N,:) - z ) ; pi(N,:) ];
  end

  % final summation (i.e. quadruple for K==2)
  pi = cumsum(a);
  z = diff(pi);
  q = ( pi(1:N-1,:) - ( pi(2:N,:) - z ) ) + ( a(2:N,:) - z );
  
  if K>0
    res = sum(q) + pi(N,:);
  else
    setround(-1)
    resinf = sum(q) + pi(N,:);
    setround(1)
    ressup = sum(q) + pi(N,:);
    res = infsup(resinf,ressup);
  end
    
  if transp
    res = res.';
  end
  
  setround(rndold)
