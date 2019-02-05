function a = permute(a,order)
%PERMUTE      Permute array dimensions
%
%Implements the Matlab function "permute", similar syntax and semantics
%

% written  09/17/08     S.M. Rump
%

  if isa(a,'intval')
    if a.complex
      a.mid = permute(a.mid,order);
      if ~isequal(a.rad,0)
        a.rad = permute(a.rad,order);
      end
    else
      a.inf = permute(a.inf,order);
      a.sup = permute(a.sup,order);
    end
  else
    error('invalid call of permute for intervals')
  end
