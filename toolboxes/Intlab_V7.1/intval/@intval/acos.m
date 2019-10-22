function y = acos(x)
%ACOS         Implements  acos(x)  for intervals
%
%   y = acos(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  pos/neg split, major revision,
%                                  improved accuracy, corrected
%                                  branchcut
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and
%                                     some improvements
% modified 09/06/07     S.M. Rump  approximate std fcts removed, exceptional arguments
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/18/08     S.M. Rump  StdFctsException ignore/NaN
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  if issparse(x)
    index = ( x==0 );
    if any(index(:))                    % treat zero indices
      INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');
      PI2 = intval(INTLAB_STDFCTS_PI.PI2INF,INTLAB_STDFCTS_PI.PI2SUP,'infsup');
      y = intval(repmat(PI2,size(x)));
      index = ~index;
      %VVVV  y(index) = acos(full(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acos(full(subsref(x,s))));
      %AAAA  Matlab bug fix
    else
      y = acos(full(x));
    end
    return
  end
      
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if x.complex
%   y = -i * log( x + sqrt(x.^2-1) );
    y = i * acosh(x);
    index = ( real(y.mid) < 0 );
    y.mid(index) = -y.mid(index);
    if rndold
      setround(rndold)
    end
    return
  end

  % input x real and full
  % real range of definition:  [-1,1]
  INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
  indexneg = ( x.inf<-1 );              % (partially) exceptional indices
  indexpos = ( x.sup>1 );               % (partially) exceptional indices
  if ~isempty(find(indexneg)) | ~isempty(find(indexpos)) % handle input out-of-range
    if INTLAB_STDFCTS_EXCPTN<=1  % out-of-range input handled as complex
      if INTLAB_STDFCTS_EXCPTN==1
        warning('ACOS: Real interval input out of range changed to be complex')
      end
      y = x;
      exceptions = indexneg | indexpos;
      %VVVV  y(exceptions) = acos(cintval(x(exceptions)));
      s.type = '()'; s.subs = {exceptions}; y = subsasgn(y,s,acos(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~exceptions;
      if any(index(:))
        %VVVV  y(index) = acos(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,acos(subsref(x,s)));
        %AAAA  Matlab bug fix
      end
      if rndold
        setround(rndold)
      end
      return
    end
    setappdata(0,'INTLAB_STDFCTS_EXCPTN_',1);
    if INTLAB_STDFCTS_EXCPTN==3    % ignore input out of range (ignore-mode)
      x.inf(indexneg) = -1;               % completely exceptional indices treated below
      x.sup(indexpos) = 1;
      Index = ( x.sup<-1 ) | ( x.inf>1 ); % completely exceptional indices
    end
  else
    Index = [];                           % make sure Index is not undefined
  end

  % input x real and full
  y = x;

  % treat non-exceptional cases
  xinf = x.inf(:);
  xsup = x.sup(:);

  % switch off warning since exceptional values may occur
  wng = warning;
  warning off

  index = ( xsup>=0 );
  if any(index)
    y.inf(index) = - asin_pos_( xsup(index) , 1 , -1 );
  end

  index = ( xsup<0 );
  if any(index)
    y.inf(index) = asin_pos_( -xsup(index) , -1 , 1 );
  end

  index = ( xinf>=0 );
  if any(index)
    y.sup(index) = - asin_pos_( xinf(index) , -1 , -1 );
  end

  index = ( xinf<0 );
  if any(index)
    y.sup(index) = asin_pos_( -xinf(index) , 1 , 1 );
  end

  if INTLAB_STDFCTS_EXCPTN==3      % ignore input out of range (ignore-mode)
    if ~isempty(find(Index))              % completely exceptional arguments to NaN
      y.inf(Index) = NaN;
      y.sup(Index) = NaN;
    end
  else                                    % any input out of range to NaN (NaN-mode)
    index = indexneg | indexpos;
    if ~isempty(find(index))              % exceptional arguments to NaN
      y.inf(index) = NaN;
      y.sup(index) = NaN;
    end
  end

  % restore warning status
  warning(wng);

  y.inf = reshape(y.inf,size(x.inf));
  y.sup = reshape(y.sup,size(x.sup));

  setround(rndold)
