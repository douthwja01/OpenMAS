function r = hull(varargin)
%HULL         interval hull
%
%   r = hull(a,b,c,...);
%
% variable parameter list
%

% only called when no parameter is interval
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  l = length(varargin);
  if l==2
    r = hull(intval(varargin{1}),varargin{2});  % calls intval function hull
  elseif l>2
    r = intval(varargin{1});
    for i=2:l
      r = hull(r,varargin{i});          % calls intval function hull
    end
  elseif l==1
    r = varargin{1};
  else
    error('hull called without parameter')
  end
  