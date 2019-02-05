function A = GregK420(n);
%GREGK420     Test matrix from Gregory/Karney
%
% Gregory/Karney test case 4.20, n>=3
%
% ( -1   2   1                      )
% (  2   0   2   1                  )
% (  1   2   0   2   1              )
% (      1   2   0   2   1          )
% (          1   2   0   2   1      )
% (              .........          )
% (              1   2   0   2   1  )
% (                  1   2   0   2  )
% (                      1   2  -1  )
%
%   A = GregK420(n);
%

% written  03/01/95     S.M. Rump
% modified 10/29/97     S.M. Rump   sparse
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if n<3
    error('dimension too small')
  end

  e=ones(n,1);
  A = spdiags( [e 2*e 0*e 2*e e], -2:2, n, n);
  A(1,1) = -1; A(n,n) = -1;
