function t = typeof(a)
%TYPEOF       Type of a
%
%   t = typeof(a)
%
%This function is to be used by typeadj for nonlinear functions. A typical
%application is
%
%   function y = f(x)
%     p = typeadj( midrad(3.14159,1e-5) , typeof(x) );
%     y = exp(p*x)-x;
%
%Then the call
%   y = f(-1.5);
%adjusts the type of p to be double by taking the midpoint, the call
%   y = f(gradient(-1.5));
%adjusts the type of p by taking the gradient of the midpoint, the call
%   y = f(gradient(intval(-1.5));
%adjusts the type of p by taking the gradient of the interval, and so forth.
%
%The table of possible values and adjustments is as follows (the first column
%gives abbreviations, the second the result of the third):
%
% d    'double'            typeof(1.4+3i);
% i    'intval'            typeof(infsup(1,2));
% g    'gradient'          typeof(gradient(4711));
% gi   'gradientintval'    typeof(gradient(midrad(2,1e-3)));
% h    'hessian'           typeof(hessian(4711));
% hi   'hessianintval'     typeof(hessian(midrad(2,1e-3)));
% s    'slope'             typeof(slopeinit(1,infsup(1,2));
% p    'polynom'           typeof(polynom([1 0 -2]));
% pi   'polynomintval'     typeof(polynom(intval([1 0 -2])));
% p    'taylor'            typeof(taylor((1:3)',5));
% pi   'taylorintval'      typeof(taylor(intval((1:3)'),5));
%
%For the conversions by typeadj, see intval\typeadj.
%

% written  10/16/98     S.M. Rump
% modified 12/06/98     S.M. Rump
% modified 12/18/02     S.M. Rump  Hessians added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/22/09     S.M. Rump  taylor added
%

% intval\typeof:  a  must be double
  t = 'double';
  