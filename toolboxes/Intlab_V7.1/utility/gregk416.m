function A = GregK416(n);
%GREGK416     Test matrix from Gregory/Karney
%
% Gregory/Karney test case 4.16, n>=3
%
% (  5  -4   1                      )
% ( -4   6  -4   1                  )
% (  1  -4   6  -4   1              )
% (      1  -4   6  -4   1          )
% (          1  -4   6  -4   1      )
% (              .........          )
% (              1  -4   6  -4   1  )
% (                  1  -4   6  -4  )
% (                      1  -4   5  )
%
%   A = GregK416(n);
%

% written  12/16/94     S.M. Rump
% modified 10/30/97     S.M. Rump   sparse
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if n<3
    error('dimension too small')
  end

  e = ones(n,1);
  A = spdiags( [e -4*e 6*e -4*e e], -2:2, n, n);
  A(1,1) = 5; A(n,n) = 5;
