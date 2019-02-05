function [res,R,p,sigma,Ms] = TransformK(p,rho,sigma,Ms)
%TRANSFORMK   Transformation res+R+sum(ps) = rho+sum(p)
%
%   [res,R,ps] = TransformK(p,rho)
%
%Input vector p may be single or double precision.
%
%Implements Algorithm 6.2 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation II: 
%    Sign, K-fold Faithful and Rounding to Nearest, Siam J. Sci. Comput., 
%    31(2):1269-1302, 2008.
%Requires (4m+3)n flops for m executions of repeat-until loop if
%  sigma,Ms not given, otherwise (4m+1)n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  if nargin==4                          % sigma and Ms already known
    kPhi = 0;
    [tau1,tau2,p,sigma,Ms] = Transform(p,rho,kPhi,sigma,Ms);
  else                                  % standard call
    [tau1,tau2,p,sigma,Ms] = Transform(p,rho);
  end
  res = tau1 + ( tau2 + sum(p) );
  R = tau2 - ( res - tau1 );
  