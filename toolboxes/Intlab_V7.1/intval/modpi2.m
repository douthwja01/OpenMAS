function [ yinf , ysup , s ] = modpi2(x)
%MODPI2       Argument reduction modulo pi/2
%
%  [ yinf , ysup , s ] = modpi2(x)
%
%For some integer k it is
%
%  x = -pi/2 + y + s*pi/4 + 2*k*pi       for even s
%  x =        -y + (s-1)*pi/4 + 2*k*pi   for odd  s
%
%with   0 <= s <= 7  and 0 <= y <= pi/4  and  yinf <= y <= ysup
%
%segment 0       1       2       3       4       5       6       7
% ÄÄÄÅÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÎÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÅÄÄ ...
%  -pi/2   -pi/4     0      pi/4    pi/2   3pi/4    pi     5pi/4   6pi/4
%
%Function modpi2 is monotone in the sense that x1<=x2, x2-x1<6 implies
%  either  yinf1 <= yinf2, s1 = s2  or  s1~=s2
%
%Works accurately thru the entire floating point range
%

% written  12/30/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  sizex = size(x);
  x = x(:);

  [ yinf , ysup , s , index ] = modpi2fast(x);

  if any(index)
    [ yinf(index) , ysup(index) , s(index) ] = modpi2large(x(index));
  end

  yinf = reshape(yinf,sizex);
  ysup = reshape(ysup,sizex);
  s = reshape(s,sizex);
    
  setround(rndold)
 


function [ zinf , zsup , s , indexpoor ] = modpi2fast(x)
% fast version of function modpi2 for not too large input
% on output, index are components with not 1 ulp result

  INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

  indexneg = ( x<0 );
  x = abs(x);
  zinf = zeros(size(x));
  zsup = zinf;
  s = zinf;
  indexpoor = ( x>2^27 );
  if all(indexpoor)
    return
  end

  k = round(x*2/pi);
  x = ( x + (-k)*INTLAB_STDFCTS_PI.PI2_1 ) + ... % exact for abs(x)<=2^27
        (-k)*INTLAB_STDFCTS_PI.PI2_2;

  setround(-1)
  zinf = x + (-k)*INTLAB_STDFCTS_PI.PI2_3SUP;
  zinf1 = ( x + (-k)*INTLAB_STDFCTS_PI.PI2_3 ) + ...
            (-k)*INTLAB_STDFCTS_PI.PI2_4SUP;
  zinf = max( zinf , zinf1 );

  setround(1)
  zsup = x + (-k)*INTLAB_STDFCTS_PI.PI2_3INF;
  zsup1 = ( x + (-k)*INTLAB_STDFCTS_PI.PI2_3 ) + ...
            (-k)*INTLAB_STDFCTS_PI.PI2_4INF;
  zsup = min( zsup , zsup1 );

  zinf1 = zinf + realmin*eps;        % rounding upwards
  indexpoor = ( zinf1<zsup ) | indexpoor;   % not 1 ulp result

  r = k - 4*floor(k/4);
  index = ( zinf>=0 );
  if any(index)                      % nonnegative output of transformation
    s(index) = mod(2*r(index)+2,8);  % zero intervals catched by indexpoor
  end
  index = ~index;
  if any(index)                      % negative output of transformation
    zinf1 = zsup(index);
    zsup(index) = -zinf(index);
    zinf(index) = -zinf1;
    s(index) = 2*r(index)+1;
  end

  if any(indexneg)                   % special treatment of negative input
    s(indexneg) = 3 - s(indexneg) + 8*(s(indexneg)>=4);
  end

  indexpoor = indexpoor | ( zinf>INTLAB_STDFCTS_PI.PI4INF );



function [ zinf , zsup , s ] = modpi2large(x)
% function modpi2 for input vector x with  1.57 <= abs(x) < inf

  INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');

  n = length(x);
  indexneg = ( x<0 );
  x = abs(x);

  % split x into mantissa and exponent; 0.5 <= f < 1
  [f,e] = log2(x);

  b = 25;                              % 2/pi given to base 2^b
  beta = 2^b;
  q = floor(e/b);
  r = e - b*q;
  F = zeros(n,3);
  for i=1:3                            % f split into 25 bits chunks
    f = beta*f;
    F(:,i) = floor(f);
    f = f - F(:,i);
  end

  % Calculate y s.t. sum(y_i*beta^(2-i)) = 2^(-r)*x*2/pi mod 4
  kmax = 6;
  if n==0
    y = conv( INTLAB_STDFCTS_PI.INV_PI2(3+q-2:3+q+kmax) , F );
  else
    y = zeros(n,kmax+5);
    Ones = ones(1,kmax+3);
    v = 3 + q*Ones + ones(n,1)*(-2:kmax);
    y(:,1:kmax+3) = (F(:,1)*Ones) .* INTLAB_STDFCTS_PI.INV_PI2(v);
    y(:,2:kmax+4) = y(:,2:kmax+4) + ...
       (F(:,2)*Ones) .* INTLAB_STDFCTS_PI.INV_PI2(v);
    y(:,3:kmax+5) = y(:,3:kmax+5) + ...
       (F(:,3)*Ones) .* INTLAB_STDFCTS_PI.INV_PI2(v);
  end
  % leading coefficient times beta^1

  % reduce y and multiply by 2^r
  y = y(:,2:end);             % leading coefficient times beta^0
  Ones = ones(1,kmax+4);
  t = floor( (2.^(r-b)*Ones) .* y );
  y = (2.^r*Ones) .* ( y - (2.^(b-r)*Ones).*t );
  y(:,1:end-1) = y(:,1:end-1) + t(:,2:end);

  % normalize mod beta, i.e. produce beta-digits
  while any( y(:)>=beta )
    R = floor( y/beta );
    y = y - beta*R;
    y(:,1:end-1) = y(:,1:end-1) + R(:,2:end);
  end

  % check result mod 4
  s = y(:,1) - 4*floor(y(:,1)/4);
  % sum(y_i*beta(-i)) <= x/pi mod 1 <= sum(y_i*beta(-i)) + 1ulp
  y = y(:,2:kmax);            % leading coefficient times beta^(-1), error 1 ulp

  % reduce y to -1/2 .. 1/2
  indexhalf = ( y(:,1)>=beta/2 );      % y >= 1/2, now 1-y stored to stay positive
  if any(indexhalf)
    y(indexhalf,:) = beta - 1 - y(indexhalf,:);     % all exact
    s(indexhalf) = s(indexhalf) + 1;
    s(s==4) = 0;
    % correction  y(:,end) = y(:,end) + 1   omitted =>
    % sum(y_i*beta(-i)) <= (1-x*2/pi) mod 1 <= sum(y_i*beta(-i)) + 1ulp
    % still valid for indices indexhalf
  end

  % multiply by pi/2
  z = conv2( y , INTLAB_STDFCTS_PI.PI2_4 );
  z = fliplr( z(:,1:kmax) .* ( ones(n,1) * ( beta.^(-2:-1:-kmax-1) ) ) );

  % calculate bounds
  setround(-1)
  zinf = sum(z,2);

  setround(1)
  yabs = ( (y(:,2)+1)/beta + y(:,1) ) / beta;
  z(:,1) = yabs*beta^(-4) + 3.15*beta^(1-kmax);
  zsup = sum(z,2);

  % correct segment
  s = 2*s + 2;
  s(indexhalf) = s(indexhalf) - 1;

  if any(indexneg)
    s(indexneg) = 11 - s(indexneg);
  end

  s(s>=8) = s(s>=8) - 8;

  setround(0)
