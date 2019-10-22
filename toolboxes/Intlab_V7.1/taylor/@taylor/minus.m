function r = minus(a,b)
%MINUS        Taylor subtraction a - b
%

% written  05/21/09     S.M. Rump
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(a,'taylor')               % non-Taylor plus Taylor
    r = b;
    if isa(a,'intval')
      r.t = intval(r.t);
    end
    r.t = -r.t;
    if isscalar(a)
      r.t(1,:) = a + r.t(1,:);
    else
      sizeb = prod(b.size);
      if prod(sizeb)==1
        r.size = size(a);
        r.t = repmat(r.t,prod(r.size));
        r.t(1,:) = a(:).' + r.t(1,:);
      else
        if ~isequal(size(a),b.size)
          error('operands of different size')
        end
        r.t(1,:) = a(:).' + r.t(1,:);
      end
    end
  elseif ~isa(b,'taylor')           % Taylor plus non-Taylor
    r = a;
    if isa(b,'intval')
      r.t = intval(r.t);
    end
    if isscalar(b)
      r.t(1,:) = r.t(1,:) - b;
    else
      sizea = prod(a.size);
      if prod(sizea)==1
        r.size = size(b);
        r.t = repmat(r.t,prod(r.size));
        r.t(1,:) = r.t(1,:) - b(:).';
      else
        if ~isequal(size(b),a.size)
          error('operands of different size')
        end
        r.t(1,:) = r.t(1,:) - b(:).';
      end
    end
  else                              % Taylor plus Taylor
    r = a;
    if isa(b.t,'intval')
      r.t = intval(r.t);
    end
    sa = prod(a.size);
    sb = prod(b.size);
    if sa==1                        % a is scalar
      if sb~=1                      % b is not scalar
        r.size = b.size;
        a.t = repmat(a.t,sb);
      end
    else                            % a is not scalar
      if sb==1                      % b is scalar
        b.t = repmat(b.t,sa);
      else
        if ~isequal(a.size,b.size)
          error('operands of different size')
        end
      end
    end
    r.t = a.t - b.t;
  end
  
  if rndold
    setround(rndold)
  end
