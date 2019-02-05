function x = circulant(x);
%CIRCULANT    Circulant matrix with first row x
%
%   A = circulant(x);
%

% written   8/27/95     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  x = toeplitz(x([1 length(x):-1:2]),x);
