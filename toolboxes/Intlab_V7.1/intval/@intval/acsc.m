function y = acsc(x)
%ACSC         Implements  acsc(x)  for intervals
%
%   y = acsc(x)
%
%interval standard function implementation
%

% written  10/16/98     S.M. Rump
% modified 12/30/98     S.M. Rump  improved speed and use asin
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy near 1
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 01/20/03     S.M. Rump  Matlab sqrt fixed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and
%                                     some improvements
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/20/08     S.M. Rump  check for zero omitted
% modified 10/20/08     S.M. Rump  StdFctsException ignore/NaN
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
    if issparse(x.mid)
      x.mid = full(x.mid);
      x.rad = full(x.rad);
    end
    y = asin( 1./x );   
    if rndold
      setround(rndold)
    end
    return
  end

  if issparse(x.inf)
    x.inf = full(x.inf);
    x.sup = full(x.sup);
  end
  % input x real and full
  % real range of definition:  [-inf,-1] and [1,inf]
  % take care for intersection( x , (-1,1) ) nonempty
  INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
  Index1 = ( abs(x.inf)<1 );           % (partially) exceptional indices
  Index2 = ( abs(x.sup)<1 ); 
  if ~isempty(find(Index1)) | ~isempty(find(Index2))   % handle input out-of-range
    if INTLAB_STDFCTS_EXCPTN<=1 % out-of-range input handled as complex
      if INTLAB_STDFCTS_EXCPTN==1
        warning('ACSC: Real interval input out of range changed to be complex')
      end
      y = x;
      index = Index1 | Index2;
      %VVVV  y(index) = asin(cintval(1./x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,asin(cintval(1./subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~index;
      if any(index(:))
        %VVVV  y(index) = asin(1./x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,asin(1./subsref(x,s)));
        %AAAA  Matlab bug fix
      end
      if rndold
        setround(rndold)
      end
      return
    end
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',1);
    if INTLAB_STDFCTS_EXCPTN==3    % ignore input out of range (ignore-mode)
      x.inf(Index1) = 1;                  % completely exceptional indices treated below
      x.sup(Index2) = -1;
      Index = ( Index1 & Index2 );        % completely exceptional indices
    end
  else
    Index = [];                           % make sure Index is not undefined
  end
  
  % input x real and full
  y = x;
  wng = warning;
  warning off
  
  % treat positive intervals
  index1 = ( x.inf>0 );
  if any(index1(:))
    y.inf(index1) = acsc_pos(x.sup(index1),-1);
    y.sup(index1) = acsc_pos(x.inf(index1),1);
  end

  % treat negative intervals
  index1 = ( x.sup<0 );
  if any(index1(:))
    y.inf(index1) = - acsc_pos(-x.sup(index1),1);
    y.sup(index1) = - acsc_pos(-x.inf(index1),-1);
  end
  
  if INTLAB_STDFCTS_EXCPTN==3      % ignore input out of range (ignore-mode)
    if ~isempty(find(Index))              % completely exceptional arguments to NaN
      y.inf(Index) = NaN;
      y.sup(Index) = NaN;
    end
  else                                    % any input out of range to NaN (NaN-mode)
    Index = Index1 | Index2;
    if ~isempty(find(Index))              % exceptional arguments to NaN
      y.inf(Index) = NaN;
      y.sup(Index) = NaN;
    end
  end

  setround(rndold)
  warning(wng)
    

  
function y = acsc_pos(x,rnd)
% local acsc for double array x>=1 with rounding corresponding to rnd
%

  y = x;

  index = ( x<1.5 );              % 1 <= x < 1.5
  if any(index(:))
    e = x(index) - 1;             % difference exact because x near 1
    setround(-rnd)
    e = sqrt_rnd( e.*(2+e) , rnd );
    setround(rnd)
    e = abs(1./e);                % take care of 1/(-0)=-inf
    y(index) = atan_pos(e(:),rnd);
  end

  index1 = ( x>1e17 );            % x > 1e17
  if any(index1(:))
    setround(-rnd)
    e = x(index1) - eps;
    setround(rnd)
    e = 1./e;
    y(index1) = atan_pos(e(:),rnd);
  end

  index = ~( index | index1 );    % 1.5 <= x <= 1e17
  if any(index(:))
    setround(-rnd)
    e = x(index);
    e = sqrt_rnd(e.*e-1,-rnd);
    setround(rnd)
    e = 1./e;
    y(index) = atan_pos(e(:),rnd);
  end
