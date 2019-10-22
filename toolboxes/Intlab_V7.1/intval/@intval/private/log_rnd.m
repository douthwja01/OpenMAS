function [x,xx] = log_rnd(x,rnd)
%LOG_RND      Rigorous bounds for log(x) according to round
%
%   y = log_rnd(x)
%
%Rounding needs not be to nearest after execution
%

% written  12/30/98     S.M. Rump
% modified 08/31/99     S.M. Rump  improved accuracy near 1, major revision
% modified 04/04/04     S.M. Rump  rounding corrected (thanks to Klaus Hildebrandt)
% modified 12/03/05     S.M. Rump  special care for x=0
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
% modified 06/01/13     S.M. Rump  two output arguments
%
  
  INTLAB_STDFCTS_LOG = getappdata(0,'INTLAB_STDFCTS_LOG');
  INTLAB_STDFCTS_LOG2 = getappdata(0,'INTLAB_STDFCTS_LOG2');

  wng = warning;
  warning off
  
  if nargout==2
    xx = x;
  end

  setround(0)
  index0 = ( x==0 );
  indexinf = ( x==inf );
  [f,e] = log2(x);
  xs = pow2( round(f.*2^14) , -14 );   % max. 14 bits of mantissa,
                                       % no bit below 2^-14

  % ensure numerical stability, avoid (0.5+eps)*2^1
  index = ( e==1 );
  if any(index(:))                      %   1 <= f < 2  for e==1
    f(index) = 2 * f(index);            % 0.5 <= f < 1  otherwise
    xs(index) = 2 * xs(index);
    e(index) = 0;
  end
  % x = f*2^e, log(x) = e*log(2) + log(xs) + log(1+(f-xs)/xs)

  logxs = log(xs);                      % rounded to nearest

  index = ( f>=xs );
  if any(index(:))

    setround(rnd)
    % positive difference
    xsindex = xs(index);
    d = ( f(index) - xsindex ) ./ xsindex;      % 0 <= d < 2^-13
    eindex = e(index);

    % rounding correct because d>=0
    % lod(1+d) = d - d^2/2 + d^3/3 -+...
    if rnd==-1                            % 0 <= err <= d^5/5 < 4.5e-17*d
      logxsindex = logxs(index);
      log1d = ((( (-d)/4 + 1/3 ).*d - 0.5 ).*d) .*d + d;
      if nargout==1
        x(index) = ( logxsindex + ( log1d + ...
                    ( (-INTLAB_STDFCTS_LOG.EPS)*abs(logxsindex) + ...
                      eindex * INTLAB_STDFCTS_LOG2.INF ) ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      else
        x(index) = logxsindex;
        xx(index) = ( log1d + ( (-INTLAB_STDFCTS_LOG.EPS)*abs(logxsindex) + ...
                      eindex * INTLAB_STDFCTS_LOG2.INF ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      end
    else                                  % 0 <= err <= d^6/6 < 4.6e-21*d
      log1d = (((( d/5 - 0.25).*d + 1/3 ).*d - 0.5 ).*d) .*d + d;
      logxsindex = logxs(index);
      if nargout==1
        x(index) = ( logxsindex + ( log1d + ...
                    ( INTLAB_STDFCTS_LOG.EPS*abs(logxsindex) + ...
                      eindex * INTLAB_STDFCTS_LOG2.SUP ) ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      else
        x(index) = logxsindex;
        xx(index) = ( log1d + ( INTLAB_STDFCTS_LOG.EPS*abs(logxsindex) + ...
                      eindex * INTLAB_STDFCTS_LOG2.SUP ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      end
    end
  end

  index = ~index;
  if any(index(:))

    setround(-rnd)
    % negative difference
    xsindex = xs(index);
    d = ( xsindex - f(index) ) ./ xsindex;      % 0 <= d < 2^-13
    eindex = e(index);

    % rounding correct because all terms nonnegative
    % lod(1-d) = - ( d + d^2/2 + d^3/3 +... )
    if rnd==-1                            % 0 <= err <= d^5/5 < 4.5e-17*d
      log1d = -( (((( 0.2002*d + 0.25).*d + 1/3 ).*d + 0.5 ).*d) .*d + d );
      setround(rnd)
      if nargout==1
        x(index) = ( logxs(index) + ( log1d + ...
                    ( (-INTLAB_STDFCTS_LOG.EPS)*abs(logxs(index)) + ...
                      eindex * INTLAB_STDFCTS_LOG2.INF ) ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      else
        x(index) = logxs(index);
        xx(index) = ( log1d + ( (-INTLAB_STDFCTS_LOG.EPS)*abs(logxs(index)) + ...
                      eindex * INTLAB_STDFCTS_LOG2.INF ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      end
    else                                  % 0 <= err <= d^5/5/(1-d)^5
      log1d = -( ((( d/4 + 1/3 ).*d + 0.5 ).*d) .*d + d );
      setround(rnd)
      if nargout==1
        x(index) = ( logxs(index) + ( log1d + ...
                    ( INTLAB_STDFCTS_LOG.EPS*abs(logxs(index)) + ...
                      eindex * INTLAB_STDFCTS_LOG2.SUP ) ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;

      else
        x(index) = logxs(index);
        xx(index) = ( log1d + ( INTLAB_STDFCTS_LOG.EPS*abs(logxs(index)) + ...
                      eindex * INTLAB_STDFCTS_LOG2.SUP ) ) + ...
                      eindex * INTLAB_STDFCTS_LOG2.APP ;
      end
    end
  end
  
  if any(index0(:))
    x(index0) = -inf;
    if nargout==2
      xx(index0) = 0;
    end
  end
  if any(indexinf(:))
    x(indexinf) = inf;
    if nargout==2
      xx(indexinf) = 0;
    end
  end

  warning(wng);
  