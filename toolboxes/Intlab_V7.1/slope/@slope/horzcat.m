function u = horzcat(varargin)
%HORZCAT      Implements  [a(1) a(2) ...]  for slopes
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  u = slope(varargin{1});

  for i=2:length(varargin)
    a = slope(varargin{i});
    u.size(2) = u.size(2) + a.size(2);
    a.size(2) = u.size(2);
    if ~isequal(a.size,u.size)
      error('dimension 1 and 3 to end must match for vertical concatenation')
    end
    u.r = [ u.r ; a.r ];   % uses that arrays are stored columnwise
    u.s = [ u.s ; a.s ];   % uses that arrays are stored columnwise
  end
