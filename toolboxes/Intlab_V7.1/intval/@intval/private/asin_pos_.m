function y = asin_pos_(x,rnd,s)
% Rigorous asin(x)+s*pi/2 for nonnegative double vector x within [0,1].
%   Rounding corresponding to rnd.
% For internal use in asin, acos

% written  08/17/99     S.M. Rump
% modified 01/20/03     S.M. Rump  Matlab sqrt fixed
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

  y = x;

  % case distinctions
  index = ( x<.75 );                   % 0 <= x < .75
  if any(index)
    xx = x(index);
    setround(-rnd)
    N = sqrt_rnd( 1 + (-xx).*xx , -rnd );
    setround(rnd)
    xrnd = xx./N;                      % to be rounded according to rnd
    switch s
      case -1                          % asin(x) - pi/2
        if rnd==-1
          y(index) = atan_pos(xrnd,rnd) - INTLAB_STDFCTS_PI.PI2SUP;
        else
          y(index) = atan_pos(xrnd,rnd) - INTLAB_STDFCTS_PI.PI2INF;
        end
      case 0                           % asin(x)
        y(index) = atan_pos(xrnd,rnd);
      case 1                           % asin(x) + pi/2
        if rnd==-1
          y(index) = atan_pos(xrnd,rnd) + INTLAB_STDFCTS_PI.PI2INF;
        else
          y(index) = atan_pos(xrnd,rnd) + INTLAB_STDFCTS_PI.PI2SUP;
        end
    end
  end

  index = ~index;
  if any(index)                        % .75 <= x <= 1
    xx = x(index);
    e = 1 - xx;                        % exact because xx close to 1
    setround(-rnd)
    x_rnd = sqrt_rnd(e+e.*xx,-rnd)./xx;     % to be rounded according to -rnd
    yrnd = atan_pos(x_rnd,-rnd);
    switch s
      case -1                          % asin(x) - pi/2
        y(index) = - yrnd;
      case 0                           % asin(x)
        setround(rnd)
        if rnd==-1
          y(index) = INTLAB_STDFCTS_PI.PI2INF - yrnd;
        else
          y(index) = INTLAB_STDFCTS_PI.PI2SUP - yrnd;
        end
      case 1                           % asin(x) + pi/2
        setround(rnd)
        if rnd==-1
          y(index) = INTLAB_STDFCTS_PI.PIINF - yrnd;
        else
          y(index) = INTLAB_STDFCTS_PI.PISUP - yrnd;
        end
    end
  end
