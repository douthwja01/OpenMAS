function r = cot(a)
%COT          Taylor cotangent  cot(a)
%
%Thanks to George Corliss for providing the Taylor expansion
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
  ct = a.t;
  N = size(a.t,2);
  r.t(1,:) = cot(a.t(1,:));
  ct(1,:) = 1+r.t(1,:).^2;     % 1+a^2
  r.t(2,:) = - ct(1,:) .* a.t(2,:);
  for i=2:K
    ct(i,:) = sum( r.t(1:i,:).*r.t(i:-1:1,:) , 1 );
    r.t(i+1,:) = - sum( ct(1:i,:).*a.t(i+1:-1:2,:).*repmat((i:-1:1)',1,N) , 1 )/i;
  end

  if rndold
    setround(rndold)
  end
