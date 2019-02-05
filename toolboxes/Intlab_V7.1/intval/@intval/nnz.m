function r = nnz(a)
%NNZ          Implements  nnz(a)  for sparse interval matrix
%
%   r = nnz(a)
%
%Functionality as in Matlab.
%

% written  08/09/02     S.M. Rump 
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if a.complex
    r = nnz(spones(a.mid)+spones(a.rad));
  else
    r = nnz(spones(a.inf)+spones(a.sup));
  end
