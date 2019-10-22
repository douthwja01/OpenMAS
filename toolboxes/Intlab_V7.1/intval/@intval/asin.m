function y = asin(x)
%ASIN         Implements  asin(x)  for intervals
%
%   y = asin(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  pos/neg split, major revision,
%                                  improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  'realstdfctsexcptnignore' added and some
%                                     improvements, tocmplx replaced by cintval
% modified 09/06/07     S.M. Rump  approximate std fcts removed, exceptional arguments
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/18/08     S.M. Rump  StdFctsException ignore/NaN
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,asin(full(sx)),m,n);
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
    y = -j * log( j*x + sqrt(1-sqr(x)) );
    if rndold
      setround(rndold)
    end
    return
  end
  
  % input x real and full
  % real range of definition:  [-1,1]
  INTLAB_STDFCTS_EXCPTN = getappdata(0,'INTLAB_STDFCTS_EXCPTN');
  indexneg = ( x.inf<-1 );               % (partially) exceptional indices
  indexpos = ( x.sup>1 );                % (partially) exceptional indices
  if ~isempty(find(indexneg)) | ~isempty(find(indexpos)) % handle input out-of-range
    if INTLAB_STDFCTS_EXCPTN<=1   % out-of-range input handled as complex
      if INTLAB_STDFCTS_EXCPTN==1
        warning('ASIN: Real interval input out of range changed to be complex')
      end
      y = x;
      index = indexneg | indexpos;
      %VVVV  y(index) = asin(cintval(x(index)));
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,asin(cintval(subsref(x,s))));
      %AAAA  Matlab bug fix
      index = ~index;
      if any(index(:))
        %VVVV  y(index) = asin(x(index));
        s.type = '()'; s.subs = {index}; y = subsasgn(y,s,asin(subsref(x,s)));
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
      index = ( x.sup<-1 ) | ( x.inf>1 ); % completely exceptional indices
    end
  else
    index = [];                           % make sure Index is not undefined
  end

  % input x real and full
  y = x;

  % in this section asin is only executed for non-exceptional cases
  
  xinf = x.inf(:);
  xsup = x.sup(:);

  % switch off warning since exceptional values may occur
  wng = warning;
  warning off

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );

  Y = asin_pos_( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 , 0 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ( xinf<0 );
  len1 = sum(IndexInfNeg);
  IndexSupPos = ( xsup>0 );

  Y = asin_pos_( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 , 0 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : end );

  if INTLAB_STDFCTS_EXCPTN==3      % ignore input out of range (ignore-mode)
    if ~isempty(find(index))              % completely exceptional arguments to NaN
      y.inf(index) = NaN;
      y.sup(index) = NaN;
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

  setround(rndold)
