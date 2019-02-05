function r = sqr(a)
%SQR          Taylor (elementwise) square  sqr(a)
%

% written  05/21/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  K1 = getappdata(0,'INTLAB_TAYLOR_ORDER') + 1;

  r = a;
  r.t(1,:) = sqr(a.t(1,:));
  for j=2:K1
    r.t(j,:) = sum(a.t(1:j,:).*a.t(j:-1:1,:),1);
  end

  if rndold
    setround(rndold)
  end
