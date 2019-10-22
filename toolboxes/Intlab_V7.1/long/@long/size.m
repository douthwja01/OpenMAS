function [m,n] = size(a,dim)
%SIZE         Implements  size(a)  for long (vectors)
%
%   [m,n] = size(a,dim)
%
% functionality as Matlab function size, second argument dim optional
%

% written  04/11/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargout==2
    [m n] = size(a.sign);
  else
    if nargin==1
      m = size(a.sign);
    else
      m = size(a.sign,dim);
    end
  end
