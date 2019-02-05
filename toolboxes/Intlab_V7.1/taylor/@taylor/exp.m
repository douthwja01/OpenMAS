function r = exp(a)
%EXP          Taylor exponential  exp(a)
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                   % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  K1 = getappdata(0,'INTLAB_TAYLOR_ORDER') + 1;
  
  r = a;
  N = size(a.t,2);
  r.t(1,:) = exp(a.t(1,:));
  for j=2:K1
    r.t(j,:) = sum( repmat((1:j-1)',1,N).*r.t(j-1:-1:1,:).*a.t(2:j,:) , 1 ) ./ (j-1);
  end

  if rndold
    setround(rndold)
  end
