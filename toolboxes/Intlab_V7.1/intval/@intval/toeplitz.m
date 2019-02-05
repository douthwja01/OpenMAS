function A = toeplitz(c,r)
%TOEPLITZ     Implements  toeplitz(c,r)  for intervals
%
%   A = toeplitz(c,r)
%
% functionality as Matlab function toeplitz
%

% written  02/20/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/07/08     S.M. Rump  help
% modified 03/12/12     S.M. Rump  real/complex mixed and radius zero
%                                    (thanks to Kolumbán Sándor, Hungary)
%

  if nargin==1
    if c.complex
      A = midrad( toeplitz(c.mid) , toeplitz(c.rad) );
    else
      A = infsup( toeplitz(c.inf) , toeplitz(c.sup) );
    end
  else
    c = intval(c);
    r = intval(r);
    if c.complex | r.complex
      c = cintval(c);
      r = cintval(r);
      if isequal(c.rad,0) & ( ~isequal(r.rad,0) )
        R = toeplitz(zeros(size(c.mid)),r.rad);
      elseif ( ~isequal(c.rad,0) ) & isequal(r.rad,0)
        R = toeplitz(c.rad,zeros(size(r.mid)));
      else
        R = toeplitz(c.rad,r.rad);
      end
      A = midrad( toeplitz(c.mid,r.mid) , R );
    else
      A = infsup( toeplitz(c.inf,r.inf) , toeplitz(c.sup,r.sup) );
    end
  end
