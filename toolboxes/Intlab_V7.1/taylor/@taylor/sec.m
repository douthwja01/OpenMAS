function r = sec(a)
%SEC          Taylor secans  sec(a)
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
 
  ct = a.t;
  st = a.t;
  r = a;
  N = size(a.t,2);
  st(1,:) = sin(a.t(1,:));
  ct(1,:) = cos(a.t(1,:));
  r.t(1,:) = 1 ./ cos(a.t(1,:));
  ct1 = -ct(1,:);
  for j=2:K+1
    at_ = a.t(2:j,:);           % some 3 % faster 
    if j~=K+1
      st(j,:) = sum( repmat((1:j-1)',1,N).*ct(j-1:-1:1,:).*at_ , 1 ) ./ (j-1);
    end
    ct(j,:) = - sum( repmat((1:j-1)',1,N).*st(j-1:-1:1,:).*at_ , 1 ) ./ (j-1);
    r.t(j,:) = sum( r.t(1:j-1,:).*ct(j:-1:2,:) , 1 ) ./ ct1;
  end

  if rndold
    setround(rndold)
  end
