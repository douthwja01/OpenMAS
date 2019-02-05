function [ C , k ] = min(A,B)
%MIN          Implements minimum of long numbers (refers only to midpoints)
%
%  C = min(A)       minimum of (column) vector A
%  [C,k] = min(A)   minimum and first index achieving minimum of A
%  C = min(A,B)     entrywise minimum of A and B, A or B may be scalar
%

% written  11/06/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
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

  if nargin==1

    I = ( A.sign<0 );
    if ~any(I)                % all nonnegative
      E = min( A.exponent );
      I = ( A.exponent==E );
      k = 0;
      while ( sum(I)>1 ) & ( k<size(A.mantissa,2) )
        k = k+1;
        E = min(A.mantissa(I,k));
        I = I & ( A.mantissa(:,k)==E );
      end
      k = find(I);
%VVVV C = A(k);
      s.type = '()'; s.subs = {k}; C = subsref(A,s);
%AAAA Matlab V5.2 bug fix
    else                      % A(index) negative
      E = max( A.exponent(I) );
      I = I & ( A.exponent==E );
      k = 0;
      while ( sum(I)>1 ) & ( k<size(A.mantissa,2) )
        k = k+1;
        E = max(A.mantissa(I,k));
        I = I & ( A.mantissa(:,k)==E );
      end
      k = find(I);
%VVVV C = A(k);
      s.type = '()'; s.subs = {k}; C = subsref(A,s);
%AAAA Matlab V5.2 bug fix
    end

  else

    C = A;
    index = ( A-B > 0 );
%VVVV C(index) = B(index);
      s.type = '()'; s.subs = {index}; C = subsasgn(C,s,subsref(B,s));
%AAAA Matlab V5.2 bug fix

  end
  
  if rndold
    setround(rndold)
  end
