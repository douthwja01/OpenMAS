function a = mag(a)
%MAG          magnitude(a)  for interval scalars, vectors and matrices
%
%   c = mag(a)
%
%On return, abs(alpha) <= c for all alpha in a
%Replaces abss. Thanks to Arnold for pointing to better name.
%

% written  01/04/05     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  if a.complex
    if isequal(a.rad,0)                 % take care of huge arrays
      a = abs(a.mid);
    else
      e = 1e-30;
      if 1+e==1-e                           % fast check for rounding to nearest
        rndold = 0;
      else
        rndold = getround;
      end
      setround(1)
      a = abs(a.mid) + a.rad;
      setround(rndold)                  % set rounding to previous value
    end
  else
    a = max( abs(a.inf) , abs(a.sup) );
  end
  