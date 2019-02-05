function q = pshift(p,beta,var)
%PSHIFT       shift of argument by beta
%
%For a polynomial p and scalar beta, scaling with repect to "var".
%
%  call:     q = pshift(p,beta,var) 
%  result:   q{x} = p{x+beta}    for x denoting the dependent variable (or var).
%
%For univariate polynomial p the parameter "var" is optional (default: dependent variable).
%
%If p does not depend on var, function is executed, possibly with warning, or an error message
%  is given depending on setting of 'AccessVariable', see polynominit
%
%Execution with error bounds in case p or beta of type intval
%

% written  07/17/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  k = size(p.e,2);                                % number of variables
  if k<=1                                         % univariate polynomial
    q = p;  
    if ( nargin==3 ) & ~isequal(var,p.v)          % polynomial not depending on specified variable
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
        case 1, warning('polynomial does not depend on specified variable')
        case 2, error('polynomial does not depend on specified variable')
      end  
      setround(rndold)
      return
    end
    p_i = p;
    if isa(p.c,'intval') | isa(beta,'intval')
      beta = intval(beta);
      q.c = intval(q.c);
      IV = 1;
    else
      IV = 0;
    end
    n = p.e;               % degree of polynomial
    q.c(n+1) = polyval(p,beta);
    for i=1:n
      p_i.e = n-i;
      if IV
        p_i.c = binom(n:-1:i,intval(i)).*p.c(1:n-i+1);      % i-th derivative of p / i!
      else
        p_i.c = binom(n:-1:i,i).*p.c(1:n-i+1);    % i-th derivative of p / i!
      end
      q.c(n+1-i) = polyval(p_i,beta);
    end
  else                                            % multivariate polynomial
    if ~ischar(var)                               % dependent variable must be string
      error('dependent variable must be string')
    end
    q = p; 
    index = find(strcmp(p.v,var));                % index of dependent variable
    if isempty(index)                             % polynomial not depending on specified variable
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
        case 1, warning('polynomial does not depend on specified variable')
        case 2, error('polynomial does not depend on specified variable')
      end  
      setround(rndold)
      return
    end
    if isa(p.c,'intval') | isa(beta,'intval')
      beta = intval(beta);
      IV = 1;
    else
      IV = 0;
    end
    n = max(p.e(:,index));                        % degree of polynomial in var
    if n==0        
      setround(rndold)
      return
    end
    index1 = 1:(index-1);
    index2 = (index+1):k;
    index12 = [index1 index2];                                       % coefficient to var^n
    q.c = q.c .* ( beta.^p.e(:,index) );
    q.e(:,index) = 0; 
    for i=1:n
      p_i = p;
      I = ( p.e(:,index)<i );                     % coefficients to disappear
      p_i.e(I,:) = [];
      p_i.c(I) = [];
      if IV
        p_i.c = p_i.c .* binom(p_i.e(:,index),intval(i)) .* ( beta.^(p_i.e(:,index)-i) );
      else
        p_i.c = p_i.c .* binom(p_i.e(:,index),i) .* ( beta.^(p_i.e(:,index)-i) );
      end
      p_i.e(:,index) = i; 
      q.e = [ q.e ; p_i.e ];
      q.c = [ q.c ; p_i.c ];
    end
    [q.e,q.c] = collect(q.e,q.c);
  end
  
  setround(rndold)
