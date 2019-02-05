function r = acosh(a)
%ACOSH        Taylor inverse cosine  acosh(a)
%

% written  06/03/09     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                   % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  K = getappdata(0,'INTLAB_TAYLOR_ORDER');

  r = a;
  dt = a.t;
  dv = dt;
  dw = dt;
  N = size(a.t,2);
  dv(1,:) = sqrt(a.t(1,:)-1);
  dw(1,:) = sqrt(a.t(1,:)+1);
  dt(1,:) = dv(1,:).*dw(1,:);         % sqrt(a-1)*sqrt(a+1)
  dv2 = 2*dv(1,:);
  dw2 = 2*dw(1,:);
  for j=2:K
    dv(j,:) = ( a.t(j,:) - sum(dv(2:j-1,:).*dv(j-1:-1:2,:),1) ) ./ dv2;
    dw(j,:) = ( a.t(j,:) - sum(dw(2:j-1,:).*dw(j-1:-1:2,:),1) ) ./ dw2;
    dt(j,:) = sum(dv(1:j,:).*dw(j:-1:1,:),1);
  end
  r.t(1,:) = acosh(a.t(1,:));
  r.t(2,:) = a.t(2,:) ./ dt(1,:);
  for j=2:K
    r.t(j+1,:) = ( j*a.t(j+1,:) - sum(repmat((1:j-1)',1,N).*r.t(2:j,:).*dt(j:-1:2,:),1) ) ./ ( j*dt(1,:) );
  end

  if rndold
    setround(rndold)
  end
