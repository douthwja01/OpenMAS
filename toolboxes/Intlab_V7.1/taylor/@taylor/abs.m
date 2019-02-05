function a = abs(a)
%ABS          Taylor absolute value abs(a)
%
%   c = abs(a)
%
%Result is convex hull of left and right derivative. This allows
%also input containing zero. The expansion
%  f(x) = f(xs) + f'(zeta)*(x-xs)
%for some zeta in ch(x,xs) is still valid
%

% written  05/21/09     S.M. Rump
%

  if ~isreal(a.t)
    a.t = abs(a.t);
    a.t(2:end,:) = repmat(NaN,size(a.t(2:end,:)));
    return
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  at = a.t(1,:);
  index = ( at<0 );
  a.t(:,index) = -a.t(:,index);
  if isa(at,'intval')
    index = in0(0,at(1,:));
    a.t(:,index) = infsup(-1,1)*a.t(:,index);
  end

  if rndold
    setround(rndold)
  end
