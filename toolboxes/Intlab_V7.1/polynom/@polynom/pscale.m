function p = pscale(p,alpha,var)
%PSCALE       scale of argument by alpha
%
%For a polynomial p and scalar alpha, scaling with repect to "var".
%
%  call:     q = pscale(p,alpha,var) 
%  result:   q{x} = p{alpha*x}    for x denoting the dependent variable (or var).
%
%For univariate polynomial p the parameter "var" is optional (default: dependent variable).
%
%If p does not depend on var, function is executed, possibly with warning, or an error message
%  is given depending on setting of 'AccessVariable', see polynominit
%
%Execution with error bounds in case p or alpha of type intval
%

% written  07/17/02     S.M. Rump
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

  if size(p.e,2)<=1                               % univariate polynomial
    if ( nargin==3 ) & ~isequal(var,p.v)          % polynomial not depending on specified variable
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
        case 1, warning('polynomial does not depend on specified variable')
        case 2, error('polynomial does not depend on specified variable')
      end
      if rndold
        setround(rndold)
      end
      return
    end
    n = p.e;                                      % degree of polynomial
    p.c = ( alpha.^(n:-1:0) ) .* p.c ;
  else                                            % multivariate polynomial
    if ~ischar(var)                               % dependent variable must be string
      error('dependent variable must be string')
    end
    index = find(strcmp(p.v,var));                % index of dependent variable
    if isempty(index)                             % polynomial not depending on specified variable
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
        case 1, warning('polynomial does not depend on specified variable')
        case 2, error('polynomial does not depend on specified variable')
      end  
      if rndold
        setround(rndold)
      end
      return
    end
    p.c = ( alpha.^p.e(:,index) ) .* p.c;
  end
  
  if rndold
    setround(rndold)
  end
