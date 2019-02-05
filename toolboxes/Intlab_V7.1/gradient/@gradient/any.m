function res = any(a,dim)
%ANY          Like Matlab function "any" for gradient
%
%Call
%
%   L = any(A)
%   L = any(A,dim)
%
%Same functionality as Matlab/any for intval quantity A
%

% written  12/06/05     S.M. Rump
%

  if nargin==1
    res = any(a.x) | any( reshape(any(a.dx,2),size(a.x)) );
  else
    res = any(a.x,dim) | any( reshape(any(a.dx,2),size(a.x)) , dim );
  end
