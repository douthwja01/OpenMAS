function y = sin_pos(x,rnd)
% rigorous sin(xs) with rounding according to rnd for 0 <= x < pi/4
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
  xs = pow2( floor(pow2(x,14)) , -14 );   % max. 14 bits of mantissa,
                                          % no bit below 2^-14
  d = x - xs;                             % 0 <= d < 2^-14
  sinxs = sin(xs);
  cosxs = cos(xs);

  % uses sin(xs+d) = sin(xs)*cos(d) + cos(xs)*sin(d)
  setround(rnd)

  if rnd==-1

    % sin(xs) >= sinxs*(1-SINEPS)
    % cos(xs) >= cosxs*(1-COSEPS)
    d2 = (-d) .* d;
    % sin(d) >= sind
    sind = d + d.*d2/6;
    % cos(d) >= 1 - d^2/2,  (-SINEPS)*cos(d) >= SINEPS*(d^2/2-1)
    y = sinxs + ...
        ( ( (sinxs.*d2)/2 + ...
            sinxs.*(INTLAB_STDFCTS_SIN.EPS.*(d.*d/2-1)) ) + ...
          cosxs.*((1-INTLAB_STDFCTS_COS.EPS).*sind) );

  else

    % sin(xs) <= sinxs*(1+SINEPS)
    % cos(xs) <= cosxs*(1+COSEPS)
    d2 = d.*d;
    % sin(d) <= sind
    sind = d + ((( d2/20 - 1 ).*d).*d).*d/6;
    % cos(d) <= 1 + cosd_
    cosd_ = (( d2/12 - 1 ).*d).*d/2;
    y = sinxs + ...
        ( ( sinxs.*cosd_ + sinxs.*INTLAB_STDFCTS_SIN.EPS.*(1+cosd_) ) + ...
            cosxs.*((1+INTLAB_STDFCTS_COS.EPS).*sind) );
  end
