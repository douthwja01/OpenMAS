function u = vertcat(varargin)
%VERTCAT      Implements  [a(1) ; a(2) ; ...]  for slopes
%

% written  12/06/98     S.M. Rump
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  u = slope(varargin{1});
  u.r = reshape(u.r,[u.size INTLAB_SLOPE.NUMVAR+1]);
  u.s = reshape(u.s,[u.size INTLAB_SLOPE.NUMVAR]);

  for i=2:length(varargin)
    a = slope(varargin{i});
    if ~isequal(a.size(2:end),u.size(2:end))
      error('dimension 2 to end must match for vertical concatenation')
    end
    u.size(1) = u.size(1) + a.size(1);
    a.r = reshape(a.r,[a.size INTLAB_SLOPE.NUMVAR+1]);
    a.s = reshape(a.s,[a.size INTLAB_SLOPE.NUMVAR]);
    u.r = [ u.r ; a.r ];
    u.s = [ u.s ; a.s ];
  end

  u.r = reshape(u.r,[prod(u.size) INTLAB_SLOPE.NUMVAR+1]);
  u.s = reshape(u.s,[prod(u.size) INTLAB_SLOPE.NUMVAR]);
