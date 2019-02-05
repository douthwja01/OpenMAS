function p = VecSum(p)
%VECSUM       Error-free vector transformation
%
%   q = VecSum(p)
%
%On return, sum(q)=sum(p) and the condition number of the sum decreases
%  by about a factor eps.
%
%Implements algorithm VecSum from
%  T. Ogita, S.M. Rump, S. Oishi: Accurate Sum and Dot Product, 
%    SIAM Journal on Scientific Computing (SISC), 26(6):1955-1988, 2005 .
%Requires 6(n-1) flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
%

  for i=2:length(p)
    [p(i),p(i-1)] = TwoSum(p(i),p(i-1));
  end

