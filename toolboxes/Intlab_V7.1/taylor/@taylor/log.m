function r = log(a)
%LOG          Taylor logarithm  log(a)
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
  r.t(1,:) = log(a.t(1,:));
  r.t(2,:) = a.t(2,:) ./ a.t(1,:);  % careful: sum(a-[])=0, not a !
  for j=3:K1
    r.t(j,:) = ( (j-1)*a.t(j,:) - sum( repmat((1:j-2)',1,N).*a.t(j-1:-1:2,:).*r.t(2:j-1,:) , 1 ) ) ...
                      ./ ((j-1)*a.t(1,:));
  end

  if rndold
    setround(rndold)
  end
