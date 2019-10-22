function p = randpoly(n,k,density,pv)
%RANDPOLY     Polynomial with random coeficients
%
%   p = randpoly(n,k,density,vars)
%
%Generates polynomial p in k variables vars{1}..vars{k} with maximum degree n
%  per variable with approximately (n+1)*k*density nonzero coefficients
%  uniformly distributed within [-1,1]. Input  vars  is a cell array of strings.
%Default of variables is 'x' for univarite and 'x1' ... 'xk' for multivariate polynomials.
%
%For example,
%  p = randpoly(5,2,.2,{'x','y'})
%produces a polynomial in variables 'x' and 'y' of maximum degree 5 per variable.
%
%Possible calls:
%
%   p = randpoly(n)                 dense univariate polynomial of degree n in 'x'
%   p = randpoly(n,var)             dense univariate polynomial of degree n in (string) var
%   p = randpoly(n,k)               polynomial in 'x1'..'xk' with density 0.5
%   p = randpoly(n,k,vars)          polynomial in vars{1}..vars{k} with density 0.5
%   p = randpoly(n,k,density)       polynomial in 'x1'..'xk' with specified density
%   p = randpoly(n,k,density,vars)  polynomial in vars{1}..vars{k} with specified density
%
%Note that  randpoly(n)  produces a univariate polynomial in 'x', whereas
%  randpoly(n,1)  produces a univariate polynomial in 'x1' and
%  randpoly(n,3)  produces a multivariate polynomial in 'x1', 'x2' and 'x3'
%To specify a univariate polynomial of degree n with prespecified density in variable 'x' use
%  randpoly(n,1,density,'x')
%

% written  09/16/00     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 02/17/08     S.M. Rump  nargin==2 corrected
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if nargin==1                                % call p = randpoly(n)
    p = polynom(random(1,n+1),'x');
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargin<=3
    if ( nargin~=2 ) | ( ~ischar(k) )
      pv{1} = 'x1';
      for i=2:k
        pv{i} = [ 'x' int2str(i) ];
      end
    end
    if nargin==2                              
      if ischar(k)                             % call p = randpoly(n,var)
        pv = k;
        k = 1;
        density = 1;
      else                                    % call p = randpoly(n,k)
        density = 0.5;
      end
    elseif iscell(density) | ischar(density)  % call p = randpoly(n,k,vars)
      pv = density;
      density = 0.5;
    end
  end
  
  L = length(pv);
  if k==1 
    if iscell(pv)
      pv = pv(1);
    else
      L = 1;
    end
  end
  if L~=k
    error('number of variables does not match number of vars')
  end
  
  pe = randint( n+1 , round((n+1)*k*density)+1 , k ) - 1;
  pc = random(size(pe,1),1);
  p = polynom(pc,pe,pv);
  
  if rndold
    setround(rndold)
  end
