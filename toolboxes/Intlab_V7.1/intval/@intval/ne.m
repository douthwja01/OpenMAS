function res = ne(a,b)
%NE           Implements  a ~= b  elementwise for intervals a and b
%
% a and b must be either both real or both complex
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump
% modified 06/24/99     S.M. Rump  multi-dimensional arrays
%                                  comparison real/complex
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if ~isa(a,'intval')
    if b.complex
      res = ( b.mid~=a ) | ( b.rad~=0 );
    else
      res = ( b.inf~=a ) | ( b.sup~=a );
    end
    return
  end

  if ~isa(b,'intval')
    if a.complex
      res = ( a.mid~=b ) | ( a.rad~=0 );
    else
      res = ( a.inf~=b ) | ( a.sup~=b );
    end
    return
  end

  if a.complex ~= b.complex
    if a.complex
      error('comparison ~= of complex and real')
    else
      error('comparison ~= of real and complex')
    end
  end

  if a.complex
    res = (a.mid~=b.mid) | (a.rad~=b.rad);
  else
    res = (a.inf~=b.inf) | (a.sup~=b.sup);
  end
