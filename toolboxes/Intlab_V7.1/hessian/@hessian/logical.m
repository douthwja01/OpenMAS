function res = logical(a)
%LOGICAL      Like Matlab function "logical" for hessian
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

  if isreal(a)                % real input
    res = logical(a.x) | reshape(any(a.dx,1),size(a.x)) | reshape(any(a.hx,1),size(a.x));
  else                        % complex input
    res = logical(real(a.x)) | logical(imag(a.x)) | ...
            reshape(any(a.dx,1),size(a.x)) | reshape(any(a.hx,1),size(a.x));
  end
  
  warning(wng)
