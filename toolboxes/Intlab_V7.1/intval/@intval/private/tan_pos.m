function y = tan_pos(x,rnd)
% rigorous tan(xs) with rounding according to rnd for 0 <= x < pi/4
% rounding remains setround(rnd) after execution
% for internal use in rigorous tan, cot

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  revision
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_TAN = getappdata(0,'INTLAB_STDFCTS_TAN');

  setround(0)
  xs = pow2( floor(pow2(x,14)) , -14 );   % max. 14 bits of mantissa,
                                          % no bit below 2^-14
  d = x - xs;                             % 0 <= d < 2^-14
  tanxs = tan(xs);                        % round to nearest

  setround(rnd)
  corr = INTLAB_STDFCTS_TAN.EPS;

  % use tan(x) = tanh(xs+d) =
  %         tanh(xs) + tan(d)*(1+tan^2(xs))/(1-tan(xs)*tan(d))
  if rnd==-1
    % tan(xs) >= tanxsinf
    tanxsinf = tanxs*(1-corr);
    % tan(d) >= tand
    tand = d + d.*d.*d/3;
    N = tanxsinf.*tand - 1;
    y = tanxs + ...
        ( tand.*( ( 1 + tanxsinf.*tanxsinf )./( -N ) ) + (-corr)*tanxs );
  else
    % tan(xs) <= tanxssup
    tanxssup = tanxs*(1+corr);
    % tan(d) <= tand,  0 <= err(tand) <= tan(IV)(d)*d^4/4! <= 9.3e-18*d
    tand = d + d .* ( d.*d/3 + 9.3e-18 );
    N = tanxssup.*tand - 1;
    y = tanxs + ...
        ( tand.*( ( 1 + tanxssup.*tanxssup )./( -N ) ) + corr*tanxs );
  end
