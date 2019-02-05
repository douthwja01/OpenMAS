function A = circulant(r)
%TOEPLITZ     Implements  circulant(r)  for intervals
%
%   A = circulant(r)
%
%Interval circulant matrix with first row r
%

% written  02/20/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if r.complex
    A = midrad( circulant(r.mid) , circulant(r.rad) );
  else
    A = infsup( circulant(r.inf) , circulant(r.sup) );
  end
