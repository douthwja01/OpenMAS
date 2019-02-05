function y = acot(x)
%ACOT         Implements  acot(x)  for intervals
%
%   y = acot(x)
%
%interval standard function implementation
%

% written  10/16/98     S.M. Rump
% modified 06/24/99     S.M. Rump  complex allowed, sparse input,
%                                  major revision, improved accuracy
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    accelaration for sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/06/07     S.M. Rump  improved performance
% modified 10/20/08     S.M. Rump  check for zero
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

  index = ( x==0 );
  y = atan(1./x);
  
  if ~isempty(find(index))                    % treat zero indices
    INTLAB_STDFCTS_PI = getappdata(0,'INTLAB_STDFCTS_PI');
    PI2 = intval(INTLAB_STDFCTS_PI.PI2INF,INTLAB_STDFCTS_PI.PI2SUP,'infsup');
    %VVVV  y(index) = PI2;
    s.type = '()'; s.subs = {index}; y = subsasgn(y,s,PI2);
    %AAAA  Matlab bug fix
  end
