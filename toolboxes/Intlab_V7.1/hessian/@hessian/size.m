function [m,n] = size(a,dim)
%SIZE         Implements  size(a)  for Hessians
%
%   [m,n] = size(a,dim)
%
% functionality as Matlab function size, second argument dim optional
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargout==2
    [m n] = size(a.x);
  else
    if nargin==1
      m = size(a.x);
    else
      m = size(a.x,dim);
    end
  end
