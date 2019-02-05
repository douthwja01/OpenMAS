function y = cosh(x)
%COSH         Implements  cosh(x)  for intervals
%
%   y = cosh(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, array input, sparse input,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/04/05     S.M. Rump  extreme values for approximate part
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/21/08     S.M. Rump  NaN input
% modified 02/19/09     S.M. Rump  sparse input
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  if issparse(x)
    index = ( x==0 );
    y = cosh(full(x));
  if ~isempty(find(index))              % treat zero indices
      %VVVV  y(index) = 1;
      s.type = '()'; s.subs = {index}; y = subsasgn(y,s,1);
      %AAAA  Matlab bug fix
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
    y = ( exp(x) + exp(-x) ) / 2;  
    if rndold
      setround(rndold)
    end
    return
  end

  if issparse(x.inf)                   % real case
    x.inf = full(x.inf);
    x.sup = full(x.sup);
  end

  y = x;

  xinf = x.inf(:);
  xsup = x.sup(:);
  indexnan = isnan(xinf) | isnan(xsup);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );
  len2 = sum(IndexSupNeg);

  Y = coshpos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.inf(IndexSupNeg) = Y( len1+1 : end );

  IndexInfNegSupPos = ~( IndexInfPos | IndexSupNeg );
  if any(IndexInfNegSupPos)
    len3 = sum(IndexInfNegSupPos);
    Y = coshpos( [  xsup(IndexInfPos) ; -xinf(IndexSupNeg) ; ...
      -xinf(IndexInfNegSupPos) ; xsup(IndexInfNegSupPos) ] , 1 );
    y.sup(IndexInfPos) = Y(1:len1);
    y.sup(IndexSupNeg) = Y( len1+1 : len1+len2 );
    y.inf(IndexInfNegSupPos) = 1;
    y.sup(IndexInfNegSupPos) = ...
      max( Y( len1+len2+1 : len1+len2+len3 ) , Y( len1+len2+len3+1 : end ) );
  else
    Y = coshpos( [ xsup(IndexInfPos) ; -xinf(IndexSupNeg) ] , 1 );
    y.sup(IndexInfPos) = Y(1:len1);
    y.sup(IndexSupNeg) = Y( len1+1 : end );
  end

  % take care of NaN input
  if ~isempty(find(indexnan))
    y.inf(indexnan) = NaN;
    y.sup(indexnan) = NaN;
  end
  setround(rndold)

  

function y = coshpos(x,rnd)
% Value cosh(x) rounded according to rnd for nonnegative double vector x

  INTLAB_STDFCTS_SINH = getappdata(0,'INTLAB_STDFCTS_SINH');

  y = x;

  % medium input
  index = ( x<709 );
  if any(index)
    [ explb , expub ] = exp_rnd(x(index));
    setround(rnd)
    if rnd==-1
      y(index) = 0.5 * ( explb + 1./expub );
    else
      y(index) = 0.5 * ( expub + 1./explb );
    end
  end

  % large input
  index = ( x>=709 );
  if any(index)
    INTLAB_STDFCTS_E = getappdata(0,'INTLAB_STDFCTS_E');
    expbnd = exp_rnd( x(index)-1 , rnd )/2;
    setround(rnd)
    if rnd==-1
      y(index) = expbnd .* INTLAB_STDFCTS_E.INF;
    else
      y(index) = expbnd .* INTLAB_STDFCTS_E.SUP;
    end
  end
