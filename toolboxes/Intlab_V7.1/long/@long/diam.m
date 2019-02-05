function C = diam(A)
%DIAM         Diameter of long (only with error term)
%
%  C = diam(A)
%

% written  11/06/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 09/20/04     S.M. Rump  output number rather than interval  
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  if ~INTLAB_LONG_ERROR
    error('diam called for long without error term')
  end

  C = sup( sup(A) - inf(A) );
