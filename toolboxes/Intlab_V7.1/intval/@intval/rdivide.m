function c = rdivide(a,b)
%RDIVIDE      Interval elementwise right division a ./ b
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  modified for infinity
% modified 06/06/98     S.M. Rump  modified for NaN+Nan*i
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    remove check for 'double'
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/03/05     S.M. Rump  sparse flag corrected
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%                                    improved performance
% modified 05/23/06     S.M. Rump  sparse Inf/NaN bug corrected in Version 7.2+
% modified 12/03/06     S.M. Rump  Sparse Bug global flag (thanks to Arnold)
% modified 10/23/07     S.M. Rump  complex numbers
% modified 02/18/09     S.M. Rump  NaN performance improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
  end
  
  % no need to take care about huge matrices: result would be almost full anyway
  % also no care necessary about previous Matlab sparse NaN bug (would be helpful, in fact)
  % make sure a full except b is scalar, b full anyway
  if isa(b,'intval')
    bcomplex = b.complex;
    if bcomplex
      b.mid = full(b.mid);
      b.rad = full(b.rad);
      makefull = ( prod(size(b.mid))~=1 );
      nanindex = isnan(b.mid) | isnan(b.rad);
    else
      b.inf = full(b.inf);
      b.sup = full(b.sup);
      makefull = ( prod(size(b.inf))~=1 );
      nanindex = isnan(b.inf) | isnan(b.sup);
    end
  else
    bcomplex = ~isreal(b);
    b = full(b);
    makefull = ( prod(size(b))~=1 );
    nanindex = isnan(b);
  end
  nanindex = sparse(nanindex);      % careful: full(True) | sparse = full
  if isa(a,'intval')
    acomplex = a.complex;
    if acomplex
      if makefull
        a.mid = full(a.mid);
        a.rad = full(a.rad);
      end
      nanindex = nanindex | isnan(a.mid) | isnan(a.rad);
    else
      if makefull
        a.inf = full(a.inf);
        a.sup = full(a.sup);
      end
      nanindex = nanindex | isnan(a.inf) | isnan(a.sup);
    end
  else
    acomplex = ~isreal(a);
    if makefull
      a = full(a);
    end
    nanindex = nanindex | isnan(a);
  end
  anynanindex = any(nanindex);
  anynanindex = any(anynanindex(:));

  ws = warning;
  warning off

  if acomplex | bcomplex                % numerator complex
    b = intval(b);                      % make sure b is interval
    if ~bcomplex                        % denominator is real
      c = a.*(1./b);
      return
    end
    x = real(b.mid);                    % denominator is complex
    y = imag(b.mid);
    setround(-1)
    Ninf = x.*x + y.*y + (-b.rad).*b.rad;
    index = ( Ninf<=0 );
    setround(1)
    Nsup = x.*x + y.*y + (-b.rad).*b.rad;
    x2 = max( x./Ninf , x./Nsup );
    y2 = max( y./Ninf , y./Nsup );
    setround(-1)
    x1 = max( x./Ninf , x./Nsup );
    y1 = max( y./Ninf , y./Nsup );
    c1 = x1 - j*y2;
    setround(1)
    c2 = x2 - j*y1;
    binv.complex = 1;
    binv.inf = [];
    binv.sup = [];
    binv.mid = c1 + 0.5*(c2-c1);
    binv.rad = abs( binv.mid - c1 ) + b.rad./Ninf;
    index = index | ( binv.rad<0 );
    if any(index(:))                          % division by zero
      binv.mid(index) = complex(NaN,NaN);
      binv.rad(index) = NaN;
    end
    binv = class(binv,'intval');
    c = a.*binv;
    if anynanindex
      c.mid(nanindex) = NaN;
      c.rad(nanindex) = NaN;                  % radius for sparse cannot be 0
    end
  else                                        % both a and b real
    c.complex = 0;
    if ~isa(a,'intval')                       % R ./ IR
      % be sure min/max works correct for zero upper bounds in b
      b.sup(b.sup==0) = -0;
      setround(-1)
      c.inf = min( a./b.inf , a./b.sup );
      setround(1)
      c.sup = max( a./b.inf , a./b.sup );
      index = ( b.inf<=0 ) & ( b.sup>=0 );    % 0 in b
      if ~isempty(find(index))
        c.inf(index) = -inf;
        c.sup(index) =  inf;
        index = index & ( a==0 );             % 0/0
        if any(index(:))
          c.inf(index) = NaN;
          c.sup(index) = NaN;
        end
        if anynanindex
          c.inf(nanindex) = NaN;
          c.sup(nanindex) = NaN;
        end
      end
    elseif ~isa(b,'intval')                     % IR ./ R
      setround(-1)
      c.inf = min( a.inf./b , a.sup./b );
      setround(1)
      c.sup = max( a.inf./b , a.sup./b );
      index = ( b==0 );                         % numerator/0
      if ~isempty(find(index))
        c.inf(index) = -inf;
        c.sup(index) =  inf;
        index = index & ( a.inf<=0 ) & ( 0<=a.sup );
        if any(index(:))
          c.inf(index) = NaN;
          c.sup(index) = NaN;
        end
        if anynanindex
          c.inf(nanindex) = NaN;
          c.sup(nanindex) = NaN;
        end
      end
    else                                        % IR ./ IR
      % be sure min/max works correct for zero upper bounds in b
      b.sup(b.sup==0) = -0;
      setround(-1)
      c.inf = min( a.inf./b.inf , a.inf./b.sup );
      c.inf = min( c.inf , a.sup./b.inf );
      c.inf = min( c.inf , a.sup./b.sup );
      setround(1)
      c.sup = max( a.inf./b.inf , a.inf./b.sup );
      c.sup = max( c.sup , a.sup./b.inf );
      c.sup = max( c.sup , a.sup./b.sup );
      index = ( b.inf<=0 ) & ( b.sup>=0 );      % 0 in b
      if ~isempty(find(index))
        if prod(size(b.inf))==1
          c.inf = -inf*ones(size(c.inf));
          c.sup = -c.inf;
        else
          c.inf(index) = -inf;
          c.sup(index) =  inf;
        end
        index = index & ( a.inf<=0 ) & ( a.sup>=0 );   % 0./0
        if any(index(:))
          c.inf(index) = NaN;
          c.sup(index) = NaN;
        end
      end
      if anynanindex
        c.inf(nanindex) = NaN;
        c.sup(nanindex) = NaN;
      end
    end
    c.mid = [];
    c.rad = [];
    c = class(c,'intval');
  end
  
  if issparse(b) & ( prod(size(c))~=1 )
    c = sparse(c);
  end

  warning(ws)
  setround(rndold)
  