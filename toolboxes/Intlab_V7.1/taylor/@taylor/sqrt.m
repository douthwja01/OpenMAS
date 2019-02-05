function r = sqrt(a)
%SQRT         Taylor square root sqrt(a)
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
  r.t(1,:) = sqrt(a.t(1,:));
  rt2 = 2*r.t(1,:);                     % almost 10 % faster
  for j=2:K1
    r.t(j,:) = ( a.t(j,:) - sum(r.t(2:j-1,:).*r.t(j-1:-1:2,:),1) ) ./ (rt2);
  end
% straight version is faster
%   for j=2:K1
%     if even(j)
%       r.t(j,:) = ( a.t(j,:)/2 - sum(r.t(2:(j/2),:).*r.t(j-1:-1:(j/2)+1,:),1) ) ./ r.t(1,:);
%     else
%       r.t(j,:) = ( a.t(j,:)/2 - sum(r.t(2:((j-1)/2),:).*r.t(j-1:-1:(j+3)/2,:),1) - ...
%                               (r.t((j+1)/2,:).^2)/2 ) ./ r.t(1,:);
%     end
%   end

  if rndold
    setround(rndold)
  end
