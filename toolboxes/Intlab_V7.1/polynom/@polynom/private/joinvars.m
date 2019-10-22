function [z,I,J] = joinvars(x,y)
%JOINVARS     joins sets of variables x and y into one set z
%
%   [z,I,J] = joinvars(x,y)
%
%On entry, string or string cell array x is a set of mutually disjoint variables,
%  and string or string cell array y is a set of mutually disjoint variables.
%On return, variables x and y are merged into z such that
%  z(I) corresponds to the set x, and z(J) corresponds to the set y.
%

% written  08/28/00     S.M. Rump
%

  % make sure, variables are stored in cell arrays (covers univariate case)
  if ~iscell(x)
    x = {x};
  end
  if ~iscell(y)
    y = {y};
  end
  
  nx = length(x);
  [z,dummy,IJ] = unique([x y]);
  I = IJ(1:nx);
  J = IJ(nx+1:end);
