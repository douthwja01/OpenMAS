function A = hankel(c,r)
%HANKEL       Implements  hankel(c,r)  for intervals
%
%   A = hankel(c,r)
%
% functionality as Matlab function toeplitz
%

% written  11/07/08     S.M. Rump
% modified 03/12/12     S.M. Rump  real/complex mixed and radius zero
%                                    (thanks to Kolumbán Sándor, Hungary)
%
% 
  if nargin==1
    if c.complex
      A = midrad( hankel(c.mid) , hankel(c.rad) );
    else
      A = infsup( hankel(c.inf) , hankel(c.sup) );
    end
  else
    c = intval(c);
    r = intval(r);
    if c.complex | r.complex
      c = cintval(c);
      r = cintval(r);
      if isequal(c.rad,0) & ( ~isequal(r.rad,0) )
        R = hankel(zeros(size(c.mid)),r.rad);
      elseif ( ~isequal(c.rad,0) ) & isequal(r.rad,0)
        R = hankel(c.rad,zeros(size(r.mid)));
      else
        R = hankel(c.rad,r.rad);
      end
      A = midrad( hankel(c.mid,r.mid) , R );
    else
      A = infsup( hankel(c.inf,r.inf) , hankel(c.sup,r.sup) );
    end
  end
