function A = GregK316(n);
%GREGK316     Example from Gregory/Karney
%
% Gregory/Karney test case 3.16, n>=2
%
% (  2  -1               )
% ( -1   2  -1           )
% (     -1   2  -1       )
% (        ........      )
% (         -1   2   1   )
% (             -1   2   )
%
%   A = GregK316(n);
%

% written   3/01/95     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if n<2
    error('dimension too small')
  end
  e = ones(n,1);
  A = spdiags( [-e 2*e -e], -1:1, n, n);
