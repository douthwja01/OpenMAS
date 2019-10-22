function  [ yinf , ysup ] = sin_(xinf,xinfinf,xinfsup,Sinf,  ...
                                 xsup,xsupinf,xsupsup,Ssup );
% rigorous sin for transformed arguments
% for internal use in sin, cos

% written  12/30/98     S.M. Rump
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

  yinf = xinf;
  ysup = xsup;

  Tinf = floor(Sinf/2);
  Tsup = floor(Ssup/2);

  % indices with result +/- 1
  setround(1)
  delta = xsup-xinf;
  index1 = ( delta >= INTLAB_STDFCTS_PI.TWOPIINF ) | ...
           ( ( Sinf==Ssup ) & ( 2*floor(Sinf/2)==Sinf) & ( xinfinf>xsupinf ) ) | ...
           ( ( Sinf==Ssup ) & ( 2*floor(Sinf/2)~=Sinf) & ( xinfinf<xsupinf ) ) | ...
           ( ( Sinf<=3 ) & ( Ssup<Sinf ) ) | ...
           ( ( Sinf>Ssup ) & ( Ssup>=4 ) );
  yinf(index1) = -1;
  ysup(index1) = 1;
  Index = ~index1;

  index = Index & ( ( Tinf<=1 ) & ( Tsup<=1 ) );  % monotonically increasing
  if any(index)

    j = index & ( Sinf==0 );
    if any(j)
      yinf(j) = - cos_pos(xinfinf(j),1);
    end
    j = index & ( Sinf==1 );
    if any(j)
      yinf(j) = - sin_pos(xinfsup(j),1);
    end
    j = index & ( Sinf==2 );
    if any(j)
      yinf(j) = sin_pos(xinfinf(j),-1);
    end
    j = index & ( Sinf==3 );
    if any(j)
      yinf(j) = cos_pos(xinfsup(j),-1);
    end

    j = index & ( Ssup==0 );
    if any(j)
      ysup(j) = - cos_pos(xsupsup(j),-1);
    end
    j = index & ( Ssup==1 );
    if any(j)
      ysup(j) = - sin_pos(xsupinf(j),-1);
    end
    j = index & ( Ssup==2 );
    if any(j)
      ysup(j) = sin_pos(xsupsup(j),1);
    end
    j = index & ( Ssup==3 );
    if any(j)
      ysup(j) = cos_pos(xsupinf(j),1);
    end

  end
  Index = Index & ~index;

  index = Index & ( ( Tinf>=2 ) & ( Tsup>=2 ) );  % monotonically decreasing
  if any(index)

    j = index & ( Ssup==4 );
    if any(j)
      yinf(j) = cos_pos(xsupsup(j),-1);
    end
    j = index & ( Ssup==5 );
    if any(j)
      yinf(j) = sin_pos(xsupinf(j),-1);
    end
    j = index & ( Ssup==6 );
    if any(j)
      yinf(j) = - sin_pos(xsupsup(j),1);
    end
    j = index & ( Ssup==7 );
    if any(j)
      yinf(j) = - cos_pos(xsupinf(j),1);
    end

    j = index & ( Sinf==4 );
    if any(j)
      ysup(j) = cos_pos(xinfinf(j),1);
    end
    j = index & ( Sinf==5 );
    if any(j)
      ysup(j) = sin_pos(xinfsup(j),1);
    end
    j = index & ( Sinf==6 );
    if any(j)
      ysup(j) = - sin_pos(xinfinf(j),-1);
    end
    j = index & ( Sinf==7 );
    if any(j)
      ysup(j) = - cos_pos(xinfsup(j),-1);
    end

  end
  Index = Index & ~index;

  index = Index & ( ( Tinf<=1 ) & ( Tsup>=2 ) );    % upper bound 1
  if any(index)

    j = index & ( Sinf==0 );
    if any(j)
      yinf(j) = - cos_pos(xinfinf(j),1);
    end
    j = index & ( Sinf==1 );
    if any(j)
      yinf(j) = - sin_pos(xinfsup(j),1);
    end
    j = index & ( Sinf==2 );
    if any(j)
      yinf(j) = sin_pos(xinfinf(j),-1);
    end
    j = index & ( Sinf==3 );
    if any(j)
      yinf(j) = cos_pos(xinfsup(j),-1);
    end

    j = index & ( Ssup==4 );
    if any(j)
      yinf(j) = min( yinf(j) , cos_pos(xsupsup(j),-1) );
    end
    j = index & ( Ssup==5 );
    if any(j)
      yinf(j) = min( yinf(j) , sin_pos(xsupinf(j),-1) );
    end
    j = index & ( Ssup==6 );
    if any(j)
      yinf(j) = min( yinf(j) , - sin_pos(xsupsup(j),1) );
    end
    j = index & ( Ssup==7 );
    if any(j)
      yinf(j) = min( yinf(j) , - cos_pos(xsupinf(j),1) );
    end

    ysup(index) = 1;

  end
  Index = Index & ~index;

  index = Index & ( ( Tinf>=2 ) & ( Tsup<=1 ) );    % lower bound -1
  if any(index)

    j = index & ( Sinf==4 );
    if any(j)
      ysup(j) = cos_pos(xinfinf(j),1);
    end
    j = index & ( Sinf==5 );
    if any(j)
      ysup(j) = sin_pos(xinfsup(j),1);
    end
    j = index & ( Sinf==6 );
    if any(j)
      ysup(j) = - sin_pos(xinfinf(j),-1);
    end
    j = index & ( Sinf==7 );
    if any(j)
      ysup(j) = - cos_pos(xinfsup(j),-1);
    end

    j = index & ( Ssup==0 );
    if any(j)
      ysup(j) = max( ysup(j) , - cos_pos(xsupsup(j),-1) );
    end
    j = index & ( Ssup==1 );
    if any(j)
      ysup(j) = max( ysup(j) , - sin_pos(xsupinf(j),-1) );
    end
    j = index & ( Ssup==2 );
    if any(j)
      ysup(j) = max( ysup(j) , sin_pos(xsupsup(j),1) );
    end
    j = index & ( Ssup==3 );
    if any(j)
      ysup(j) = max( ysup(j) , cos_pos(xsupinf(j),1) );
    end

    yinf(index) = -1;

  end
