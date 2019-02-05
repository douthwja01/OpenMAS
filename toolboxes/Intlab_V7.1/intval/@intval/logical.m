function res = logical(a)
%LOGICAL      Like Matlab function "logical" for intval
%
%Call
%
%   L = logical(A)
%
%Same functionality as Matlab/logical for intval quantity A
%

% written  12/06/05     S.M. Rump
%

  wng = warning;
  warning off

  if a.complex                % complex input
    res = logical(real(a.mid)) | logical(imag(a.mid)) | logical(a.rad);
  else                        % real input
    res = logical(a.inf) | logical(a.sup);
  end
 
  warning(wng)
