function [y,ysup] = erf_rnd(x,rnd)
% input x real non-negative column vector
% rnd  -1  y = lower bound for erf(x)
%       1  y = upper bound for erf(x)
%      []  [y,ysup] inclusion of erf(x)
% rounding may be altered after leaving erf_rnd
%

% written  05/30/13     S.M. Rump
%

  x2 = hex2num('4017744f8f74e94b');     % ~5.86: erf(x)>pred(1) for x>=x2
  
  y = x;
  if isempty(rnd)
    ysup = x;
  end
  
  index = ( x<=0.5 );                   % first method
  if any(index)                         % x in [0,0.5]
    if isempty(rnd)
      [y(index),ysup(index)] = erf_rnd1(x(index),rnd,6);
    else
      y(index) = erf_rnd1(x(index),rnd,6);
    end
  end
  Index = index;                        % store finished indices
  
  index = ( ~Index ) & ( x<x2 );        % second method
  if any(index)                         % x in (0.5,x2)
    y_index = 1 - erfc_rnd2(x(index));
    if rnd==-1
      y(index) = y_index.inf;
    elseif rnd==1
      y(index) = y_index.sup;
    else
      y(index) = y_index.inf;
      ysup(index) = y_index.sup;
    end
  end
  
  index = ( x>=x2 );
  if any(index)                         % x in [x2,inf]
    if isempty(rnd)                     % inclusion [y,ysup] of erf(x)
      y(index) = hex2num('3fefffffffffffff');   % pred(1)
      ysup(index) = 1;
    elseif rnd==1                       % y upper bound for erf(x)
      y(index) = 1;
    else                                % y lower bound for erf(x)
      y(index) = hex2num('3fefffffffffffff');   % pred(1)
    end
  end
  