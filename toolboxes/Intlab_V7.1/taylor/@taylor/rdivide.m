function r = rdivide(a,b)
%RDIVIDE      Taylor elementwise right division a ./ b (same as a/b)
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

  if ~isa(a,'taylor')           % non-taylor / taylor
    r = b;
    if isa(a,'intval')
      r.t = intval(r.t);
    end
    m = prod(size(a));
    n = prod(b.size);
    if m==1                     % non-taylor scalar / taylor
      r.t(1,:) = a  ./ b.t(1,:);
      for j=2:K1
        r.t(j,:) = ( - sum(b.t(2:j,:).*r.t(j-1:-1:1,:),1) ) ./ b.t(1,:);
      end
    else                        % non-taylor array / taylor
      n = prod(b.size);
      if n==1                   % non-taylor array / taylor scalar
        r.t = repmat(r.t,1,m);
        r.t(1,:) = a ./ b.t(1);
        for j=2:K1
          r.t(j,:) = ( - sum(repmat(b.t(2:j),1,m).*r.t(j-1:-1:1,:),1) ) ./ b.t(1,:);
        end
      else                      % non-taylor array / taylor array
        if ~isequal(size(a),b.size)
          error('Taylor division : dimensions not compatible')
        end
        r.t(1,:) = a  ./ b.t(1,:);
        for j=2:K1
          r.t(j,:) = ( - sum(b.t(2:j,:).*r.t(j-1:-1:1,:),1) ) ./ b.t(1,:);
        end
      end
    end
  elseif ~isa(b,'taylor')       % taylor / non-taylor
    r = a;
    if isa(b,'intval')
      r.t = intval(r.t);
    end
    m = prod(size(b));
    if m==1                     % taylor scalar / non-taylor
      r.t = r.t/b;
    else                        % taylor array / non-taylor
      n = prod(a.size);
      if n==1                   % taylor array / non-taylor scalar
        r.t = repmat(a.t,1,m) ./ repmat(b(:).',K1,1);
      else                      % taylor array / non-taylor array
        if ~isequal(size(b),a.size)
          error('Taylor division : dimensions not compatible')
        end
        r.t = a.t ./ repmat(b(:).',K1,1);
      end
    end
  else                          % both factors taylor
    m = prod(a.size);
    n = prod(b.size);
    if m==1                     % taylor scalar ./ taylor
      r = b;
      if isa(a,'intval')
        r.t = intval(r.t);
      end
      for j=1:K1
        r.t(j,:) = ( a.t(j) - sum(b.t(2:j,:).*r.t(j-1:-1:1,:),1) ) ./ b.t(1,:);
      end
    else                        % taylor array ./ taylor
      r = a;
      if isa(b.t,'intval')
        r.t = intval(r.t);
      end
      if n==1                   % taylor array ./ taylor scalar
        for j=1:K1
          r.t(j,:) = ( a.t(j,:) - sum(repmat(b.t(2:j),1,m).*r.t(j-1:-1:1,:),1) ) ./ b.t(1);
        end
      else                      % taylor array ./ taylor array
        if ~isequal(a.size,b.size)
          error('Taylor multiplication : dimensions not compatible')
        end
        for j=1:K1
          r.t(j,:) = ( a.t(j,:) - sum(b.t(2:j,:).*r.t(j-1:-1:1,:),1) ) ./ b.t(1,:);
        end
      end
    end
  end
  
  if rndold
    setround(rndold)
  end
