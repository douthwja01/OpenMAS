function I = find(p)
%FIND         Index set of nonzero coefficients of univariate polynomial
%
%   I = find(p)
%
%For example P=polynom([1 0 -2 1]), I=find(P) results in
%  polynom P[x] = 
%      1.0000  x^3  
%     -2.0000  x    
%      1.0000       
%  I =
%       3     1     0
%

% written  10/04/02     S.M. Rump 
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if length(p.v)>1
    error('find only for univariate polynomials')
  end
  I = p.e-find(p.c)+1;
