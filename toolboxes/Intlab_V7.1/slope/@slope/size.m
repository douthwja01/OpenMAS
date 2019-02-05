function [m,n] = size(a,dim)
%SIZE         Implements  size(a)  for slope
%
%   [m,n] = size(a,dim)
%
% functionality as Matlab function size, second argument dim optional
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargout==2
    [m,n] = size(zeros(a.size));
  elseif nargin==1
    m = a.size;
  else
    m = size(zeros(a.size),dim);
  end
