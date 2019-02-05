function p = ptrans(p,a,b,c,d,var)
%PTRANS       transformation of argument
%
%For a polynomial p and scalars alpha,beta,a,b,c,d, transformation of variable "var".
%
%  q = ptrans(p,alpha,beta,var)  yields   q{x} = p{alpha*x+beta}
%  q = ptrans(p,a,b,c,d,var)     yields   q{c} = p{a} and q{d} = p{b}   transformation of [a,b] to [c,d]
%
%For univariate polynomial p the parameter "var" is optional (default: dependent variable).
%
%Execution with error bounds in case one of p,alpha,beta,a,b,c or d of type intval
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

  checkarg = 1;
  if nargin<=4                                  % call ptrans(p,a,b) or ptrans(p,a,b,var)
    if nargin==3       
      if size(p.e,2)==1                         % univariate polynomial
        c = p.v;
        checkarg = 0;
      else
        error('invalid call: variable not specified.')
      end
    end                                       
    if checkarg                                 % check variable var (is c)
      index = find(strcmp(c,p.v));
      if isempty(index)                         % polynomial not depending on specified variable
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
    end
    p = pscale(pshift(p,b,c),a,c);
  else
    if nargin==5       
      if size(p.e,2)==1                         % univariate polynomial
        var = p.v;
        checkarg = 0;
      else
        error('invalid call: variable not specified.')
      end
    end                                       
    if checkarg
      index = find(strcmp(var,p.v));
      if isempty(index)                         % polynomial not depending on specified variable
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
    end
    if isa(p.c(1),'intval') | isa(a,'intval') | isa(b,'intval') | isa(c,'intval') | isa(d,'intval')
      a = intval(a);
      c = intval(c);
    end
    alpha = (b-a)/(d-c);
    beta = (a*d-b*c)/(d-c);
    p = pscale(pshift(p,beta,var),alpha,var);
  end
  
  if rndold
    setround(rndold)
  end
