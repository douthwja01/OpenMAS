function y = atan(x)
%ATAN         Implements  atan(x)  for intervals
%
%   y = atan(x)
%
%interval standard function implementation
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,atan(full(sx)),m,n);
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
    y = log( 2./(1-j*x) - 1) / (2*j);
    if rndold
      setround(rndold)
    end
    return
  end

  y = x;

  xinf = x.inf(:);
  xsup = x.sup(:);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );

  Y = atan_pos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ( xinf<0 );
  len1 = sum(IndexInfNeg);
  IndexSupPos = ( xsup>0 );

  Y = atan_pos( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : end );

  setround(rndold)
  