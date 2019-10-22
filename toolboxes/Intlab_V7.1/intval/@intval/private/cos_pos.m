function y = cos_pos(x,rnd)
% rigorous cos(xs) with rounding according to rnd for 0 <= x < pi/4
% rounding remains setround(rnd) after execution
% for internal use in rigorous sin, cos

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  revision
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_SIN = getappdata(0,'INTLAB_STDFCTS_SIN');
  INTLAB_STDFCTS_COS = getappdata(0,'INTLAB_STDFCTS_COS');

  setround(0)
  xs = pow2( floor(pow2(x,14)) , -14 );  % max. 14 bits of mantissa,
                                         % no bit below 2^-14
  d = x - xs;                            % 0 <= d < 2^-14
  sinxs = sin(xs);
  cosxs = cos(xs);

  % uses cos(xs+d) = cos(xs)*cos(d) - sin(xs)*sin(d)

  if rnd==-1                         % result rounded downwards

    % cos(xs) >= cosxs*(1-COSEPS)
    % sin(xs) <= sinxssup
    setround(1)                      % rounding upwards
    sinxssup = sinxs*(1+INTLAB_STDFCTS_SIN.EPS);
    % sin(d) <= sind
    sind = d + (((( d.*d/20 - 1 ).*d).*d).*d)/6;
    % cos(d) >= 1 - d^2/2

    setround(-1)                     % rounding downwards
    y = cosxs + ...
        ( ( cosxs.*((-d).*d)/2 + ...
            cosxs.*( INTLAB_STDFCTS_COS.EPS.*(d.*d/2-1) ) ) + ...
          (-sinxssup).*sind );

  else                               % result rounded upwards

    % cos(xs) <= cosxs*(1+COSEPS)
    % sin(xs) >= sinxsinf
    setround(-1)                     % rounding downwards
    sinxsinf = sinxs*(1-INTLAB_STDFCTS_SIN.EPS);
    % sin(d) >= sind
    sind = d + (((-d).*d).*d)/6;
    % cos(d) <= 1 + cosd_
    setround(1)                      % rounding upwards
    cosd_ = ( d.*d/12 - 1 ) .* d.*d/2;

    y = cosxs + ...
        ( ( cosxs.*cosd_ + ...
            cosxs.*( INTLAB_STDFCTS_COS.EPS.*(1+cosd_) ) ) + ...
          (-sinxsinf).*sind );

  end
