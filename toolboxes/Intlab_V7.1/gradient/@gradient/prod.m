function c = prod(a,dim)
%PROD         Implements  prod(a,dim)  for gradients
%
%   c = prod(a,dim)
%
% functionality as Matlab function prod for matrices, parameter dim optional
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  [m n] = size(a);
  if nargin==1,
    if m==1
      dim=2;
    else
      dim=1;
    end
  end

  if dim==1
    c = ones(1,n);
    for i=1:m
%VVVV c = c .* a(i,:);
      s.type = '()'; s.subs = {i,':'}; c = c .* subsref(a,s);
%AAAA Matlab V5.2 bug fix
    end
  else
    c = ones(m,1);
    for i=1:n
%VVVV c = c .* a(:,i);
      s.type = '()'; s.subs = {':',i}; c = c .* subsref(a,s);
%AAAA Matlab V5.2 bug fix
    end
  end
  
  if rndold
    setround(rndold)
  end
