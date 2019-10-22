function [p,q] = bandwidth(A)
%BANDWIDTH    Upper and lower bandwidth of matrix A
%
%   A   = 0  for  i-j > p  or j-i > q
%    ij
%
%    [p,q] = bandwidth(A)
%

% written  10/16/98     S.M. Rump
% modified 09/21/02     S.M. Rump  one output argument corrected
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if A.complex
    [p1,q1] = bandwidth(A.mid);
    [p2,q2] = bandwidth(A.rad);
  else
    [p1,q1] = bandwidth(A.inf);
    [p2,q2] = bandwidth(A.sup);
  end

  p = max(p1,p2);
  q = max(q1,q2);
  
  if nargout<=1
    p = max(abs(p),abs(q));
  end
