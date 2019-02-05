function [y,ysup] = erfc_rnd(x,rnd)
% input x real non-negative column vector
% rnd  -1  y = lower bound for erf(x)
%       1  y = upper bound for erf(x)
%      []  [y,ysup] inclusion of erf(x)
% rounding may be altered after leaving erfc_rnd

% written  05/30/13     S.M. Rump
%

  xmax = hex2num('403b369a6244e684');   % ~27.2: erfc(x)<subrealmin fuer x>=xmax
  
  y = x;
  if isempty(rnd)
    ysup = x;
  end
  
  index = ( x<-0.5 );                   % Use 2-erfc(-x)
  if any(index)                         % x in [-inf,-0.5)
    if isempty(rnd)
      [yindex,ysupindex] = erfc_rnd(-x(index),-rnd);
      setround(-1)
      y(index) = 2 - ysupindex;
      setround(1)
      ysup(index) = 2 - yindex;
    else
      yindex = erfc_rnd(-x(index),-rnd);
      setround(rnd)
      y(index) = 2 - yindex;
    end
  end
  Index = index;                        % store finished indices
  
  index = ( ~Index ) & ( x<0 );         % first method
  if any(index)                         % x in [-0.5,0)
    if isempty(rnd)
      [yindex,ysupindex] = erf_rnd1(-x(index),rnd,6);
      setround(-1)
      y(index) = 1 + yindex;
      setround(1)
      ysup(index) = 1 + ysupindex;
    else
      yindex = erf_rnd1(-x(index),rnd,6);
      setround(rnd)
      y(index) = 1 + yindex;
    end
  end
  Index = Index | index;                % store finished indices
  
  index = ( ~Index ) & ( x<=0.5 );      % first method
  if any(index)                         % x in [0,0.5]
    if isempty(rnd)
      [yindex,ysupindex] = erf_rnd1(x(index),rnd,6);
      setround(-1)
      y(index) = 1 - ysupindex;
      setround(1)
      ysup(index) = 1 - yindex;
    else
      yindex = erf_rnd1(x(index),-rnd,6);
      setround(rnd)
      y(index) = 1 - yindex;
    end
  end
  Index = Index | index;                % store finished indices
  
  index = ( ~Index ) & ( x<=7 );        % second method
  if any(index)                         % x in (0.5,7]
    yindex = erfc_rnd2(x(index));
    if rnd==-1
      y(index) = yindex.inf;
    elseif rnd==1
      y(index) = yindex.sup;
    else
      y(index) = yindex.inf;
      ysup(index) = yindex.sup;
    end
  end
  Index = Index | index;                % store finished indices
  
  index = ( ~Index ) & ( x<=10 );       % third method
  if any(index)                         % x in (7,10]
    if isempty(rnd)
      [y(index),ysup(index)] = erfc_rnd3(x(index),rnd,16);
    else
      y(index) = erfc_rnd3(x(index),rnd,16);
    end
  end
  Index = Index | index;                % store finished indices
  
  index = ( ~Index ) & ( x<=xmax );     % third method
  if any(index)                         % x in (10,xmax]
    if isempty(rnd)
      [y(index),ysup(index)] = erfc_rnd3(x(index),rnd,10);
    else
      y(index) = erfc_rnd3(x(index),rnd,10);
    end
  end
  
  index = ( x>xmax );
  if any(index)                         % x in (xmax,inf]
    if isempty(rnd)                     % inclusion [y,ysup] of erfc(x)
      y(index) = 0;
      ysup(index) = hex2num('0000000000000001');    % subrealmin
    elseif rnd==1                       % y upper bound for erfc(x)
      y(index) = hex2num('0000000000000001');    % subrealmin
    else                                % y lower bound for erfc(x)
      y(index) = 0;
    end
  end
  