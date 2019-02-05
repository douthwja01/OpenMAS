function res = relerr(x,y);
%RELERR       Relative error between two numbers, vectors or matrices
% Below eps ( = 1e-15 )  switch to absolute error
%
%   res = relerr(x,y);
%
%Alternative call   res = relerr(X)   for interval quantity X is
%  equivalent to    res = relerr(X.inf,X.sup)
%

% written   6/12/95     S.M. Rump
% modified 11/11/97     S.M. Rump  x or y may be scalar
% modified 11/30/97     S.M. Rump  x may be NaN
% modified 03/01/01     S.M. Rump  interval input
% modified 11/16/01     S.M. Rump  denominator abs(x)+abs(y)
% modified 09/22/02     S.M. Rump  constant 1e-20 replaced by 1e-15
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    input converted to full to avoid Matlab sparse bugs
% modified 11/24/02     S.M. Rump  res=0 for one argument
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/08/05     S.M. Rump  NaN corrected
% modified 01/03/07     S.M. Rump  sparse input, redesign
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 11/13/10     S.M. Rump  result for inf or NaN
%

  if nargin==1
    if issparse(x)
      [m,n] = size(x);
      res = sparse([],[],[],m,n);
    else
      res = zeros(size(x));
    end
    res(isinf(x) | isnan(x)) = inf;     % make sure relerr of inf or NaN is infinity
    return
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  eps = 1e-15;
  if issparse(x) | issparse(y)
    [m,n] = size(x);
    N = abs(x)+abs(y);
    index = find(N);
    res = N;
    if any(index(:))
      res(index) = abs(x(index)-y(index))./N(index);
    end
    [i,j,s] = find(x);
    small = sparse(i,j,abs(s)<eps,m,n);
    [m,n] = size(y);
    [i,j,s] = find(y);
    small = small | sparse(i,j,abs(s)<eps,m,n);
    diff = abs(x-y);
    [m,n] = size(diff);
    [i,j,s] = find(diff);
    small = small | sparse(i,j,abs(s)<eps,m,n);
    res(small) = diff(small);
  else
    diff = abs(x-y);
    bool = (diff<eps) | (abs(x)<eps) | (abs(y)<eps) ;
    res = double(bool);
    if any(bool(:))
      b0=find(bool==0); b1=find(bool==1);
      scalarx = prod(size(x))==1;
      scalary = prod(size(y))==1;
      res(b1) = diff(b1);
      if scalarx
        res(b0) = abs( (x    -y(b0)) ./ (abs(x)    +abs(y(b0))) );
      elseif scalary
        res(b0) = abs( (x(b0)-y    ) ./ (abs(x(b0))+abs(y)    ) );
      else
        res(b0) = abs( (x(b0)-y(b0)) ./ (abs(x(b0))+abs(y(b0))) );
      end;
    else
      res = abs( (x-y) ./ (abs(x)+abs(y)) );
    end
  end
  
  index = find( isnan(x) | isnan(y) );
  if any(index(:))
    res(index) = NaN;
  end
  
  if rndold
    setround(rndold)
  end
