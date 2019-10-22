function y = erfc_rnd2(x)
% input x real non-negative column vector
% where  x0~.958: [x0,succ(x0)] contains root of erf^(V)
%        x1~2.02: [x1,succ(x1)] contains root of erf^(V)
% erf^(IV) monotonically increasing in (x0,x1], i.e. erf^(IV)(xi)<=erf^(IV)(x1)
% rounding may be altered after leaving erf_rnd2
%

% written  05/30/13     S.M. Rump
%

  setround(0)
  
  % factorLB <= 2/sqrt(pi) <= factorUB
  INTLAB_STDFCTS_ERF = getappdata(0,'INTLAB_STDFCTS_ERF');
  factorLB = INTLAB_STDFCTS_ERF.TWO_SQRTPIINF;
  factorUB = INTLAB_STDFCTS_ERF.TWO_SQRTPISUP;  % ~ 1.12
  
  eps_erfc = 1.085e-16;                     % erfc(x) in erfc_data(nu)*(1+/-eps)
  erfc_data = getappdata(0,'INTLAB_INTVAL_ERFCDATA');   % erfc_data loaded
 
  phi = 1024;
  xs = floor(phi*x)/phi;                    % x = xs + delta
  delta = x - xs;                           % 0 <= delta < x1/phi
  erfc_xs = erfc_data((xs-0.5)*phi+1);      % fl-pt approximation
  
  % correction xs -> xs+delta
  X = infsup(xs,x).^2;
  Delta = intval(delta);
  f13 = exp(-X);
  f2 = 4*exp(-(xs+Delta/2).^2);
  corr = Delta/6.*(f13.inf+f2+f13.sup) - Delta.^5/2880 .* f13 .* ((16*X-48).*X+12);
  y = erfc_xs + ( midrad(0,eps_erfc)*erfc_xs - infsup(factorLB,factorUB)*corr );
