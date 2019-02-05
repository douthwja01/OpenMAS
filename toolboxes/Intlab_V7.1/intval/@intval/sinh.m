function y = sinh(x)
%SINH         Implements  sinh(x)  for intervals
%
%   y = sinh(x)
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
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/04/05     S.M. Rump  extreme values for approximate part
% modified 09/06/07     S.M. Rump  approximate std fcts removed
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  if issparse(x)
    [ix,jx,sx] = find(x);
    [m,n] = size(x);
    y = sparse(ix,jx,sinh(full(sx)),m,n);
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
    y = ( exp(x) - exp(-x) ) / 2;  
    setround(rndold)
    return
  end

  y = x;

  xinf = x.inf(:);
  xsup = x.sup(:);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );
  len2 = sum(IndexSupNeg);

  Y = sinhpos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ~IndexInfPos;
  len1 = sum(IndexInfNeg);
  IndexSupPos = ~IndexSupNeg;
  len2 = sum(IndexSupPos);

  Y = sinhpos( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : len1+len2 );
  
  setround(rndold)


function y = sinhpos(x,rnd)
% rigorous sinh(x) for nonnegative double vector x with
% rounding corresponding to rnd

  INTLAB_STDFCTS_SINH = getappdata(0,'INTLAB_STDFCTS_SINH');

  y = x;

  % small input
  index = ( x<8 );
  if any(index)
    setround(0)
    xx = x(index);
    [f,e] = log2(xx);
    bits = 14 + min(e,0);
    xs = pow2( floor(f.*2.^bits) , e-bits );  % max. 14 bits of mantissa,
                                              % no bit below 2^-14
    d = xx - xs;                              % 0 <= d < 2^-14*x < 2^-11
    sinhxs = sinh(xs);                        % round to nearest

    % bounds for exp(xs)
    INTLAB_STDFCTS_EXP = getappdata(0,'INTLAB_STDFCTS_EXP');
    setround(-1)
    expinf = exp(xs)*(1-INTLAB_STDFCTS_EXP.EPS);
    setround(1)
    expsup = exp(xs)*(1+INTLAB_STDFCTS_EXP.EPS);
    expinf(xs==0) = 1;                    % gives accuracy for small x
    expsup(xs==0) = 1;

    % use sinh(xx) = sinh(xs+d) = sinh(xs)*cosh(d) + cosh(xs)*sinh(d)
    setround(rnd)
    corr = INTLAB_STDFCTS_SINH.EPS;
    if rnd==-1
      % coshxs <= cosh(xs)
      coshxs = 0.5*( expinf + 1./expsup );  % rounded downwards
      dd = d.*d;
      % 0 <= err(sinhd) <= sinh(d)*d^4/4! < 6e-8*d^3/6
      % 0 <= err(coshd) <= sinh(d)*d^5/5! < 4.8e-8*d^4/24
      y(index) = sinhxs + ...
          ( ( sinhxs .* ( -corr + (1-corr)*dd/2.*( 1 + dd/12 ) ) + ...
              coshxs.*d.*dd/6 ) + ...
            coshxs.*d ...
          );
    else
      % coshxs >= cosh(xs)
      coshxs = 0.5*( expsup + 1./expinf );  % correctly rounded according to rnd
      dd = d.*d;
      y(index) = sinhxs + ...
          ( ( sinhxs .* ( corr + (1+corr)*dd/2.*( 1 + dd/12.*( 1+4.8e-8 ) ) ) + ...
              coshxs.*d.*dd/6.*( 1+6e-8 ) ) + ...
            coshxs.*d ...
          );
    end
  end

  % medium input
  index = ( ~index ) & ( x<709 );         % 8 <= x < 709
  if any(index)
    xx = x(index);
    exprnd = exp_rnd(xx,rnd);
    setround(rnd)
    y(index) = 0.5 * ( exprnd + (-1)./exprnd );
  end

  % large input
  index = ( x>=709 );
  if any(index)
    INTLAB_STDFCTS_E = getappdata(0,'INTLAB_STDFCTS_E');
    exprnd = exp_rnd( x(index)-1 , rnd )/2;
    setround(rnd)
    if rnd==-1
      y(index) = exprnd .* INTLAB_STDFCTS_E.INF;
    else
      y(index) = exprnd .* INTLAB_STDFCTS_E.SUP;
    end
  end
