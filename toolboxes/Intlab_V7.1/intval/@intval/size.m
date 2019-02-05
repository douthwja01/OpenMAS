function [m,n] = size(a,dim)
%SIZE         Implements  size(a)  for intervals
%
%   [m,n] = size(a,dim)
%
% functionality as Matlab function size, second argument dim optional
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if a.complex
    if nargout==2
      [m n] = size(a.mid);
    else
      if nargin==1
        m = size(a.mid);
      else
        m = size(a.mid,dim);
      end
    end
  else
    if nargout==2
      [m n] = size(a.inf);
    else
      if nargin==1
        m = size(a.inf);
      else
        m = size(a.inf,dim);
      end
    end
  end
