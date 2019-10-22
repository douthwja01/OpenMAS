function res = in0(a,b)
%IN0          Implements  a in int(b)  entrywise for intervals a and b
%
%  res = in0(a,b)
%
%Result is scalar, vector or matrix depending on a and b.
%For components equal to 1,  corresponding components of  a  are
%  definitely included in the interior of the corresponding component of b
%For b being vector or matrix, "a" may be scalar
%

% written  10/16/98     S.M. Rump
% modified 11/29/98     S.M. Rump  result matrix
% modified 09/02/00     S.M. Rump  rounding unchanged after use, complex case corrected
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    remove check for 'double'
%                                    take care of huge arrays
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/04/05     S.M. Rump  tocmplx replaced by cintval
% modified 03/07/10     S.M. Rump  size check
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if ~isa(b,'intval')
    error('attempt to check  in(a,b)  for b no interval')
  end

  if ( prod(size(a))==1 ) & ( prod(size(b))~=1 )
    % a = a.*ones(size(b));
    a = repmat(full(a),size(b));
  end

  if ~isa(a,'intval')
    a = intval(a);
  end
  if ~isequal(size(a),size(b))
    error('attempt to check  in(a,b)  for a and b of different size')
  end

  if b.complex            % b is complex
    if ~a.complex
      a = cintval(a);
    end
    setround(-1)
    R1 = real(b.mid) - real(a.mid);
    I1 = imag(b.mid) - imag(a.mid);
    setround(1)
    R2 = real(b.mid) - real(a.mid);
    I2 = imag(b.mid) - imag(a.mid);
    R = max( abs(R1) , abs(R2) );
    I = max( abs(I1) , abs(I2) );
    d = abs( R + sqrt(-1)*I );
    if isequal(a.rad,0)
      res = ( d < b.rad );
    else
      res = ( (d+a.rad) < b.rad );
    end
  else                    % b is real
    if a.complex
      error('attempt to check  in(a,b)  for complex a, real b')
    end
                          % both a and b real
    res = ( b.inf < a.inf ) & ( a.sup < b.sup );
  end

  setround(rndold)
  