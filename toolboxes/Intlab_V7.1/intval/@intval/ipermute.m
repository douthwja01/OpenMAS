function a = ipermute(a,order)
%IPERMUTE     Inverse permute array dimensions
%
%Implements the Matlab function "ipermute", similar syntax and semantics
%

% written  09/17/08     S.M. Rump
%

  if isa(a,'intval')
    if a.complex
      a.mid = ipermute(a.mid,order);
      if ~isequal(a.rad,0)
        a.rad = ipermute(a.rad,order);
      end
    else
      a.inf = ipermute(a.inf,order);
      a.sup = ipermute(a.sup,order);
    end
  else
    error('invalid call of permute for intervals')
  end
