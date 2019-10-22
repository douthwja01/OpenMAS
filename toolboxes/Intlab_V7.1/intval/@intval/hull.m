function r = hull(varargin)
%HULL         Interval hull
%
%   r = hull(a,b,c,...);
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
% modified 02/20/09     S.M. Rump  radius 0 and hull(NaN,*) = NaN
% modified 05/28/09     S.M. Rump  typo
% modified 08/17/12     S.M. Rump  treatment of NaN
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  l = length(varargin);
  if l==2
    
    a = varargin{1};
    b = varargin{2};
    if isempty(a)
      r = b;  
      setround(rndold)
      return
    end
    if isempty(b)
      r = a;        
      setround(rndold)
      return
    end
    if prod(size(a))>1             % check whether a or b is scalar
      if prod(size(b))==1          % b is scalar
        if issparse(a)
          b = b*spones(a);         % preserve sparsity
        else
          b = b*ones(size(a));
        end
      end
    else                           % a is scalar
      if prod(size(b))>1
        if issparse(b)
          a = a*spones(b);         % preserve sparsity
        else
          a = a*ones(size(b));
        end
      end
    end
    
    if ~isequal(size(a),size(b))
      error('interval hull called with non-matching dimensions')
    end
    
    if ~isa(a,'intval')
      a = intval(a);
    end
    if ~isa(b,'intval')
      b = intval(b);
    end
    
    if a.complex | b.complex     % a or b is complex
      if ~a.complex
        a = cintval(a);
      end
      if ~b.complex
        b = cintval(b);
      end
      r.complex = 1;
      r.inf = [];
      r.sup = [];
      if issparse(a)
        a.rad = sparse(a.rad);
      end
      if issparse(b)
        b.rad = sparse(b.rad);
      end
      
      d = b.mid - a.mid;
      setround(1)
      dist = abs( abs(real(d)) + sqrt(-1)*abs(imag(d)) );
      r.mid = a.mid + (0.5*(dist+b.rad-a.rad)./dist).*d;
      r.rad = eps*abs(dist) + a.rad + abs(r.mid-a.mid)*(1+eps);
      
      % entire interval a in b for some indices
      index = ( dist+a.rad <= b.rad );
      if any(index(:))
        r.mid(index) = b.mid(index);
        if isequal(b.rad,0)
          r.rad(index) = 0;
        else
          r.rad(index) = b.rad(index);
        end
      end
      
      % entire interval b in a for some indices
      index = ( dist+b.rad <= a.rad );
      if any(index(:))
        r.mid(index) = a.mid(index);
        if isequal(a.rad,0)
          r.rad(index) = 0;
        else
          r.rad(index) = a.rad(index);
        end
      end
      
    else                         % both a and b are real
      r.complex = 0;
      r.inf = min(a.inf,b.inf);
      r.sup = max(a.sup,b.sup);
      if issparse(r.inf)
        [i,j] = find(isnan(a.inf) | isnan(a.sup) | isnan(b.inf) | isnan(b.sup));
        [m n] = size(a.inf);
        r.inf = r.inf + sparse(i,j,NaN,m,n);
        r.sup = r.sup + sparse(i,j,NaN,m,n);
      else
        index = isnan(a.inf) | isnan(a.sup) | isnan(b.inf) | isnan(b.sup);
        if any(index(:))
          r.inf(index) = NaN;
          r.sup(index) = NaN;
        end
      end
      r.mid = [];
      r.rad = [];
    end
    
    r = class(r,'intval');
    
  elseif l>2
    r = intval(varargin{1});
    for i=2:l
      r = hull(r,varargin{i});          % calls intval function hull
    end
  elseif l==1
    r = varargin{1};
  else
    error('hull called without parameter')
  end
        
  setround(rndold)
