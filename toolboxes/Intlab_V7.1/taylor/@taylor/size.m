function [m,n] = size(a,dim)
%SIZE         Implements  size(a)  for Taylor
%
%   [m,n] = size(a,dim)
%
% functionality as Matlab function size, second argument dim optional
%

% written  05/21/09     S.M. Rump
%

  if nargout==2
    m = a.size(1);
    n = a.size(2);
  else
    if nargin==1
      m = a.size;
    else
      m = a.size(dim);
    end
  end
