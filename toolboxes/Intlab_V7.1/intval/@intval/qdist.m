function res = qdist(a,b)
%QDIST        Implements  q(a,b)  metrical distance
%  name  qdist  is used to avoid ambiguities with variable  q
%
%     res = qdist(a,b)          max( abs(inf(a)-inf(b)) , abs(sup(a)-sup(b)) )
%
% for real intervals      max( abs(inf(a)-inf(b)) , abs(sup(a)-sup(b)) )
% for complex intervals   qdist(real(a),real(b)) + qdist(imag(a),imag(b))
%

% written  10/16/98     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(a,'intval')
    a = intval(a);
  end
  if ~isa(b,'intval')
    b = intval(b);
  end

  if ~a.complex & ~b.complex
    res = max( abs(a.inf-b.inf) , abs(a.sup-b.sup) ) ;
  else
    res = qdist(real(a),real(b)) + qdist(imag(a),imag(b)) ;
  end
  
  if rndold
    setround(rndold)
  end
