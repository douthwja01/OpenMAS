function [tau,p] = ExtractVector(p,sigma)
%EXTRACTVECTOR  Extract vector into higher and lower order part w.r.t. sigma
%
%   [tau,ps] = ExtractVector(p,sigma)
%
%On return, tau + sum(ps_i) = sum(p_i). Input vector p may be
%single or double precision.
%
%Implements Algorithm 3.3 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, SIAM J. Sci. Comput., 31(1):189-224, 2008.
%Requires 4n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  q = ( sigma + p ) - sigma;
  p = p - q;
  tau = sum(q);
  