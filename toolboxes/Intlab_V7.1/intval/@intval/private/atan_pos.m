function y = atan_pos(x,rnd)
% Rigorous atan(x) for nonnegative double vector x, rounding corresponding
%   to rnd
% For internal use in atan, asin, acos, asec, acsc

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  extreme input, improved speed and accuracy
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

  y = x;

  % transformation of large arguments
  index = ( x>=4 );
  if any(index(:))
    setround(-rnd)
    xx = x(index);
    xx = (-2) ./ ( 1./xx - xx );       % with correct rounding -rnd
    yy = atan_pos_small(xx,-rnd);
    setround(rnd)
    if rnd==-1
      y(index) = INTLAB_STDFCTS_PI.PI2INF - yy/2;
    else
      y(index) = INTLAB_STDFCTS_PI.PI2SUP - yy/2;
    end
  end

  index = ~index;
  if any(index(:))
    y(index) = atan_pos_small(x(index),rnd);
  end




function y = atan_pos_small(x,rnd)
% Rigorous atan(x) for nonnegative double vector x with  0 <= x < 4
% rounding corresponding to rnd

  INTLAB_STDFCTS_ATAN = getappdata(0,'INTLAB_STDFCTS_ATAN');

  setround(0)                               % maximum first 15 bits, no bit
  xs = pow2( floor(pow2(x,13)) , -13 );     %   below 2^-13
  d = x - xs;                               % 0 <= d < 2^-15*x
  atanxs = atan(xs);

  % atan(xs+d) = atan(xs) + atan(d/(1+x*xs))
  setround(-rnd)
  E = 1 + x.*xs;
  setround(rnd)
  E = d./E;                                 % 0 <= E < d/(1+x*xs) < 2^-13

  if rnd==-1                         % 0 <= err <= E^5/5 <= d^4/5*d < 5e-17*d
    % atanE <= atan(E)
    atanE = E + (((-E).*E).*E)/3;
    y = atanxs + ( atanE + (-INTLAB_STDFCTS_ATAN.EPS)*atanxs );
 else
    % atanE >= atan(E)
    atanE = ((( E.*E/5 + (-1)/3 ).*E).*E).*E + E;
    y = atanxs + ( atanE + INTLAB_STDFCTS_ATAN.EPS*atanxs );
  end
