function p = pderiv(p,k,var)
%PDERIV       k-th derivative of polynomial w.r.t. to dependent variable 'var'
%
%For univariate polynomials,
%
%   ps = pderiv(p)        first derivative p' of p
%   ps = pderiv(p,k)      k-th derivative of p
%
%First derivative of univariate p also implemented by ctranspose ( p' )
%
%
%For multivariate polynomials,
%
%   ps = pderiv(p,var)    first derivative w.r.t. dependent variable 'var'
%   ps = pderiv(p,var,k)  k-th derivative w.r.t. dependent variable 'var'
%

% written  09/27/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if size(p.e,2)<=1        % univariate polynomial
    n = p.e;               % degree of polynomial
    if nargin==1, k=1; end
    if k==0      
      if rndold
        setround(rndold)
      end
      return
    end
    if k==1
      if n==0              % constant polynomial
        p.c = typeadj( 0 , typeof(p.c) );
      else                 % degree at least 1
        p.c = ( n:-1:1 ) .* p.c(1:n);
      end
    elseif k>n
      p.c = typeadj( 0 , typeof(p.c) );
    else
      f = -ones(k,n-k+1);
      f(1,:) = n:-1:k;
      p.c = prod(cumsum(f)) .* p.c(1:n-k+1);
    end;
    p.e = length(p.c)-1;
  else                     % multivariate polynomial
    if ~ischar(k)          % dependent variable must be string
      error('dependent variable must be string')
    end
    index = find(strcmp(p.v,k));   % index of dependent variable
    if isempty(index)      % polynomial not depending on specified variable
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
      case 1, warning('polynomial does not depend on specified variable')
      case 2, error('polynomial does not depend on specified variable')
      end
      p.e = zeros(1,size(p.e,2));
      p.c = typeadj( 0 , typeof(p.c) );  
      if rndold
        setround(rndold)
      end
      return
    end
    if nargin==2
      k = 1;
    else
      k = var;
    end
    if k==0        
      if rndold
        setround(rndold)
      end
      return
    end
    I = ( p.e(:,index)<k );    % coefficients to disappear
    p.e(I,:) = [];
    if isempty(p.e)            % derivative is zero
      p.e = zeros(1,size(p.e,2));
      p.c = typeadj( 0 , typeof(p.c) );  
      if rndold
        setround(rndold)
      end
      return
    end
    p.c(I) = [];
    f = -ones(size(p.e,1),k);
    f(:,1) = p.e(:,index);
    p.c = prod(cumsum(f,2),2) .* p.c;
    p.e(:,index) = p.e(:,index)-k;
  end
  
  if rndold
    setround(rndold)
  end
