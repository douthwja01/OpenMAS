function y = tanh(x)
%TANH         Implements  tanh(x)  for intervals
%
%   y = tanh(x)
%
%interval standard function implementation
%

% written  10/16/98     S.M. Rump
% modified 06/24/99     S.M. Rump  complex allowed, sparse input
% modified 08/31/99     S.M. Rump  improved accuracy for small arguments,
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
    y = sparse(ix,jx,tanh(full(sx)),m,n);
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
    y = sinh(x)./cosh(x);  
    setround(rndold)
    return
  end

  y = x;

  xinf = x.inf(:);
  xsup = x.sup(:);

  IndexInfPos = ( xinf>=0 );
  len1 = sum(IndexInfPos);
  IndexSupNeg = ( xsup<=0 );

  Y = tanh_pos( [ xinf(IndexInfPos) ; -xsup(IndexSupNeg) ] , -1 );
  y.inf(IndexInfPos) = Y(1:len1);
  y.sup(IndexSupNeg) = -Y( len1+1 : end );

  IndexInfNeg = ( xinf<0 );
  len1 = sum(IndexInfNeg);
  IndexSupPos = ( xsup>0 );

  Y = tanh_pos( [ -xinf(IndexInfNeg) ; xsup(IndexSupPos) ] , 1 );
  y.inf(IndexInfNeg) = -Y(1:len1);
  y.sup(IndexSupPos) = Y( len1+1 : end );
  
  setround(rndold)


function y = tanh_pos(x,rnd)
% rigorous tanh(x) for nonnegative double vector x with
% rounding corresponding to rnd

  INTLAB_STDFCTS_TANH = getappdata(0,'INTLAB_STDFCTS_TANH');

  y = x;

  % small input
  index = ( x<1 );
  if any(index)
    setround(0)
    xx = x(index);
    xs = pow2( floor(pow2(xx,14)) , -14 );   % max. 14 bits of mantissa,
                                             % no bit below 2^-14
    d = xx - xs;                             % 0 <= d < 2^-14
    tanhxs = tanh(xs);

    % use tanh(xx) = tanh(xs+d) =
    %         tanh(xs) + tanh(d)*(1-tanh^2(xs))/(1+tanh(xs)*tanh(d))
    corr = INTLAB_STDFCTS_TANH.EPS;
    % d-d^3/3 <= tanh(d) <= d-d^3/3+err,  0 <= err <= 2/15*d^5 < 1.9e-18*d

    if rnd==-1
      setround(1)
      % tanhxssup >= tanh(xs)
      tanhxssup = tanhxs*(1+corr);
      % tanhd >= tan(d)
      tanhdsup = d + d.*( ((-d).*d)/3 + 1.9e-18 );
      N = 1 + tanhxssup.*tanhdsup;
      setround(-1)
      tanhdinf = d + (((-d).*d).*d)/3;
      % numerator >= 1 - tanh^2(1) > .4
      y(index) = tanhxs + ...
          ( tanhdinf.*( ( 1 + (-tanhxssup).*tanhxssup )./N ) + (-corr)*tanhxs );
    else
      setround(-1)
      % tanhxsinf <= tanh(xs)
      tanhxsinf = tanhxs*(1-corr);
      tanhdinf = d + (((-d).*d).*d)/3;
      N = 1 + tanhxsinf.*tanhdinf;
      setround(1)
      tanhdsup = d + d.*( ((-d).*d)/3 + 1.9e-18 );
      y(index) = tanhxs + ...
          ( tanhdsup.*( ( 1 + (-tanhxsinf).*tanhxsinf )./N ) + corr*tanhxs );
    end
  end
  index = ( ~index );

  % large input                  % x >= 1
  if any(index)
    expx = exp_rnd( 2*x(index) , rnd );
    setround(rnd)
    y(index) = 1 + (-2) ./ ( expx + 1 );
  end
