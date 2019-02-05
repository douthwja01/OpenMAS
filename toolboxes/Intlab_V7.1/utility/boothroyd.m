function [A,Ainv] = Boothroyd(n)
%ZIELKE       Ill-conditioned matrix with checkerboard sign pattern
%
% Boothroyd matrix, possibly with inverse
%
%   res = boothroyd(n)   or   [A,Ainv] = boothroyd(n)
%
%Exactly representable until n=20. The 2-norm condition numbers are
%
%     n    cond_2
% ----------------
%     1  1.00e+000
%     2  1.79e+001
%     3  4.88e+002
%     4  1.82e+004
%     5  7.86e+005
%     6  3.69e+007
%     7  1.83e+009
%     8  9.38e+010
%     9  4.96e+012
%    10  2.68e+014
%    11  1.47e+016
%    12  8.20e+017
%    13  4.63e+019
%    14  2.63e+021
%    15  1.51e+023
%    16  8.75e+024
%    17  5.10e+026
%    18  2.98e+028
%    19  1.76e+030
%    20  1.04e+032
%

% written   09/28/94     S.M. Rump
% modified  02/29/08     S.M. Rump  Comment on maximum n added
%

  A = zeros(n);
  for i=1:n
    for j=1:n
      A(i,j) = binom(n+i-1,i-1)*n*binom(n-1,n-j)/(i+j-1);
    end
  end
  if nargout==2
    Ainv = A;
    for i=1:n
      for j=1:n
        if rem(i+j,2)==1, Ainv(i,j)=-Ainv(i,j); end
      end
    end
  end

