function Y = fletcher(X)
%R. Fletcher: Distillation column test problem,
%   in More (ed.): A Collection of Nonlinear Model Problems,
%   Lectures in Applied Mathematics, pp. 723-762, 26(1990).
%
%   y = fletcher(x)
%
%Handles three different problems:
%     1  Hydrocarbon-6 problem, 29 unknowns
%     2  Hydrocarbon-20 problem, 99 unknowns
%     3  Methanol-8 problem, 31 unknowns
%
%The call
%   xs = fletcher('init1')
%generates the approximation given in the above paper for the first problem,
% accordingly parameter 'init2' and 'init3' for the second and third problem,
% respectively.
%Fletcher's functions are evaluated by
%   y = fletcher(x)
%with the function choosen by the dimension of x. A verified inclusion
% of a solution can be calculated by
%   X = verifynlss('fletcher',xs)
%

% written  09/23/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ischar(X)              % return initial values
    Y = initialvalues(X);  
    if rndold
      setround(rndold)
    end
    return
  end
  
  switch length(X)
    case 29,  index = 1;    % Hydrocarbon-6 problem, 29 unknowns
    case 99,  index = 2;    % Hydrocarbon-20 problem, 99 unknowns
    case 31,  index = 3;    % Methanol-8 problem, 31 unknowns
    otherwise
      error('invalid input to fletcher')
  end
  
  
  % get input parameters
  [ n,m,                    ... % input dimensions,
      k,                    ... % stage of distillation
      a,b,c,                ... % Antoine constants
      alpha,alphap,alphapp, ... % liquid enthalpy constants
      beta,betap,betapp,    ... % vapor anthalpy constants
      fl,fv,                ... % constants to maintain material balance
      tf,bb,d,q,            ... % constants to maintain heat balance
      pi_                   ... % stage pressures
  ] = initdata(index);
  
  % input variables (take care of starting index zero)
  x_ = reshape(X(1:n*m),m,n)';  % x(i,j) = x_(i+1,j)
  t_ = X(n*m+1:n*m+n);          % t(i) = t_(i+1)
  v_ = X(n*m+n+1:n*m+2*n-1);    % v(i) = v_(i+1)
  
  % auxiliary quantities
  l = typeadj(zeros(n-1,1),typeof(X));
  l(1:k) = v_(1:k) + bb;
  l(k+1:n-1) = v_(k+1:n-1) - d;
  
  y_ = typeadj(zeros(n,m),typeof(X));   % y(i,j) = y_(i+1,j)
  for j=1:m
    y_(:,j) = exp(a(j)+b(j)./(c(j)+t_))./pi_.*x_(:,j);
  end
  
  h_ = typeadj(zeros(n,1),typeof(X));   % h(i) = h_(i+1)
  for i=1:n
    h_(i) = sum( x_(i,:)'.*( ( alphapp*t_(i) + alphap ) *t_(i) + alpha ) );
  end
  
  H_ = typeadj(zeros(n,1),typeof(X));   % H(i) = H_(i+1)
  for i=1:n
    H_(i) = sum( y_(i,:)'.*( ( betapp*t_(i) + betap ) *t_(i) + beta ) );
  end
  
  % heat content of the feed
  hf = sum( fl .* ( ( alphapp*tf + alphap ) *tf + alpha ) );
  Hf = sum( fv .* ( ( betapp*tf + betap ) *tf + beta ) );
  
  % set of equations
  Y = X;
  
  I = 0;
  J = 1:m;
  Y(I+J) = l(1)*x_(2,:) - v_(1)*y_(1,:) - bb*x_(1,:);
  I = I+m;
  for i=1:n-2
    Y(I+J) = v_(i)*y_(i,:) + l(i+1)*x_(i+2,:) - v_(i+1)*y_(i+1,:) - l(i)*x_(i+1,:);
    if i==k
      Y(I+J) = Y(I+J) + fl;
    end
    if i==k+1
      Y(I+J) = Y(I+J) + fv;
    end
    I = I+m;
  end
  Y = .01*Y;
  Y(I+J) = y_(n-1,:) - x_(n,:);
  I = I+m;
  
  Y(I+(1:n)) = sum(y_,2) - 1;
  I = I+n;
  
  Y(I+1) = 1e-6 * ( l(1)*h_(2) + q - v_(1)*H_(1) - bb*h_(1) );
  I = I+1;
  J = 1:n-2;
  Y(I+J) = 1e-6 * ( v_(J).*H_(J) + l(J+1).*h_(J+2) - v_(J+1).*H_(J+1) - l(J).*h_(J+1) );
  Y(I+k) = Y(I+k) + 1e-6*hf;
  Y(I+k+1) = Y(I+k+1) + 1e-6*Hf ;
  
  % end of function fletcher
    
  if rndold
    setround(rndold)
  end

  
  
  % input data and constants to problem(index)
  function ...
    [ n,m,                    ... % input dimensions
      k,                    ... % stage of distillation
      a,b,c,                ... % Antoine constants
      alpha,alphap,alphapp, ... % liquid enthalpy constants
      beta,betap,betapp,    ... % vapor anthalpy constants
      fl,fv,                ... % constants to maintain material balance
      tf,bb,d,q,            ... % constants to maintain heat balance
      pi_                   ... % stage pressures
  ] = initdata(index)
  
  switch index
    
    case 1                        % Hydrocarbon-6 problem, 29 unknowns
      n = 6;
      m = 3;
      k = 2;
      a = [ 9.647 ; 9.953 ; 9.466 ];
      b = [ -2998 ; -3448.1 ; -3347.25 ];
      c = [ 230.66 ; 235.88 ; 215.31 ];
      alpha = [ 0 ; 0 ; 0 ];
      alphap = [ 37.6 ; 48.2 ; 45.4 ];
      alphapp = [ 0 ; 0 ; 0 ];
      beta = [ 8425 ; 9395 ; 10466 ];
      betap = [ 24.2 ; 35.6 ; 31.9 ];
      betapp = [ 0 ; 0 ; 0];
      fl = [ 30 ; 30 ; 40 ];
      fv = [ 0 ; 0 ; 0 ];
      tf = 100;
      bb = 40;
      d = 60;
      q = 2500000;
      pi_ = ones(6,1);
      
    case 2                        % Hydrocarbon-20 problem, 99 unknowns
      n = 20;
      m = 3;
      k = 9;
      a = [ 9.647 ; 9.953 ; 9.466 ];
      b = [ -2998 ; -3448.1 ; -3347.25 ];
      c = [ 230.66 ; 235.88 ; 215.31 ];
      alpha = [ 0 ; 0 ; 0];
      alphap = [ 37.6 ; 48.2 ; 45.4 ];
      alphapp = [ 0 ; 0 ; 0 ];
      beta = [ 8425 ; 9395 ; 10466 ];
      betap = [ 24.2 ; 35.6 ; 31.9 ];
      betapp = [ 0 ; 0 ; 0];
      fl = [ 30 ; 30 ; 40 ];
      fv = [ 0 ; 0 ; 0 ];
      tf = 100;
      bb = 40;
      d = 60;
      q = 2500000;
      pi_ = ones(20,1);
      
    case 3                        % Methanol-8 problem, 31 unknowns
      n = 8;
      m = 2;
      k = 2;
      a = [ 18.5751 ; 18.3443 ];
      b = [ -3632.649 ; -3841.2203 ];
      c = [ 239.2 ; 228 ];
      alpha = [ 0 ; 0 ];
      alphap = [ 15.97 ; 18.1 ];
      alphapp = [ .0422; 0 ];
      beta = [ 9566.67 ; 10834.67 ];
      betap = [ -1.59; 8.74 ];
      betapp = [ .0422 ; 0];
      fl = [ 451.25 ; 684.25 ];
      fv = [ 0 ; 0 ];
      tf = 89;
      bb = 693.37;
      d = 442.13;
      q = 8386200;
      pi_ = [ 1210 ; 1200 ; 1190 ; 1180 ; 1170 ; 1160 ; 1150 ; 1140 ];
      
    otherwise
      error('invalid index in fletcher initdata')
  end
  
  
  function x = initialvalues(index)
  % initial values to problem(index)
  
  switch index
    
    case 'init1'                  % Hydrocarbon-6 problem, 29 unknowns
      x = [ 0 .2 .9 0 .2 .8 .05 .3 .8 ...
          .1 .3 .6 .3 .5 .3 .6 .6 0 ...
          100 100 100 100 100 100 ...
          300 300 300 300 300]';
      
    case 'init2'                  % Hydrocarbon-20 problem, 99 unknowns
      x = [ 0 .3 1 0 .3 .9 .01 .3 .9 ...
          .02 .4 .8 .05 .4 .8 .07 .45 .8 ...
          .09 .5 .7 .1 .5 .7 .15 .5 .6 ...
          .2 .5 .6 .25 .6 .5 .3 .6 .5 ...
          .35 .6 .5 .4 .6 .4 .4 .7 .4 ...
          .42 .7 .3 .45 .75 .3 .45 .75 .2 ...
          .5 .8 .1 .5 .8 0 ...
          100*ones(1,20) 300*ones(1,19) ]';
      
    case 'init3'                  % Methanol-8 problem, 31 unknowns
      x = [ .09203 .908 .1819 .8181 .284 .716 .3051 .6949...
        .3566 .6434 .468 .532 .6579 .3421 .8763 .1237 ...
          120 110 100 88 86 84 80 76 ...
          886.37 910.01 922.52 926.45 935.56 952.83 975.73 ]';
      
    otherwise
      error('invalid index in fletcher initialvalues')
  end
