function [ y , ysup ] = exp_rnd(x,rnd)
%EXP_RND      Rigorous bounds for exp(x) and input vector x according to rnd
%
%   [ y , ysup ] = exp_rnd(x,rnd)
%
%If specified with two output arguments, lower and upper bound is computed
%independent of rnd, otherwise y is rounded according to rnd.
%Rounding need not be to nearest after leaving exp_rnd
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  extreme input, improved accuracy,
%                                  major revision
% modified 02/14/01     S.M. Rump  improved accuracy
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_EXP = getappdata(0,'INTLAB_STDFCTS_EXP');

  infsup = ( nargout==2 );
  if infsup
    rnd = -1;
  end

  setround(0)
  xint = fix(x);                            % integer part of x, chopped
  xfrac = x - xint;                         % fractional part, -1 < xfrac < 1

  indexlarge = ( xint>709 );                % care for exceptions
  indexsmall = ( xint<-744 );
  index = indexlarge | indexsmall | isnan(x);
  xint(index) = 0;
  xfrac(index) = 0;

  xs = pow2( floor(pow2(xfrac,14)) , -14 );  % max. 14 bits of mantissa of xfrac,
                                             % no bit below 2^-14
  d = xfrac - xs;                            % 0 <= d < 2^-14, exactly repr.
  expxs = exp(xs);
  factor = INTLAB_STDFCTS_EXP.EPS;
  Exint = INTLAB_STDFCTS_EXP.POW(745+xint);
  Exinteps = INTLAB_STDFCTS_EXP.POWSUP(745+xint);

  % general case, exp(xfrac) = exp(xs)*exp(d)
  if infsup          % 0 <= err <= exp(d)*d^4/4! < 0.2501*d*d^3/3!

    % calculate upper bound
    setround(1)
    % exp(d)*expxs <= corr + expxs
    corr = (( ( 1+0.2501*d ) .* expxs.*d/3 + expxs ).*d/2 + expxs ).*d;

    % exp(xs) <= expxs*(1+EPS)
    % exp(xfrac) = exp(d)*exp(xs) <= (corr+expxs)*(1+EPS)
    ysup = ( corr + corr*factor ) + expxs*factor ;

    % exp(xfrac) <= expxs+ysup,  exp(x) = exp(xfrac)*exp(xint)
    ysup = expxs.*Exint + ( (expxs+ysup).*Exinteps + ysup.*Exint );
    ysup(x==0) = 1;

  end

  setround(rnd)
  if rnd==-1
    Exinteps = INTLAB_STDFCTS_EXP.POWINF(745+xint);
    corr1 = 1;
  else
    corr1 = 1 + 0.2501*d;
  end

  % exp(d)*expxs ~ corr + expxs    subject to rounding
  corr = (( corr1 .* expxs.*d/3 + expxs ).*d/2 + expxs ).*d;

  % exp(xs)  in  expxs*(1+/-EPS)
  % exp(xfrac) = exp(d)*exp(xs)  in  (corr+expxs)*(1+rnd*EPS)
  y = ( corr + (corr*rnd)*factor ) + (expxs*rnd)*factor;

  % exp(xfrac) ~ expxs+y,  exp(x) = exp(xfrac)*exp(xint)
  y = expxs.*Exint + ( (expxs+y).*Exinteps + y.*Exint );
  y(x==0) = 1;


  % large or small input, exceptions
  if any(index)

    if any(indexlarge)
      if infsup
        y(indexlarge) = realmax;
        ysup(indexlarge) = inf;
      else
        if rnd==-1
          y(indexlarge) = realmax;
        else
          y(indexlarge) = inf;
        end
      end
    end

    if any(indexsmall)
      INTLAB_INTVAL_ETA = realmin*eps;  % smallest positive denormalized fl-pt
      if infsup
        y(indexsmall) = 0;
        ysup(indexsmall) = INTLAB_INTVAL_ETA;
      else
        if rnd==-1
          y(indexsmall) = 0;
        else
          y(indexsmall) = INTLAB_INTVAL_ETA;
        end
      end
    end

    index = isnan(x);
    if any(index)
      y(index) = NaN;
      if infsup
        ysup(index) = NaN;
      end
    end

  end

  index = ( xint<-720 );
  if any(index)
    y(index) = max( y(index) , 0 );
  end
