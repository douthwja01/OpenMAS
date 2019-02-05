function res = lssresidual(A,x,b)
%LSSRESIDUAL    Improved inclusion of residual b-A*x or I-R*A
%Typical calls are
%
%   res = lssresidual(A,x,b)
%
%Input   A     nxn real or complex point matrix (dense or sparse)
%        x     some approximation to A\b
%        b     real or complex point or interval vector or nxk matrix
%Output  res   approximation of b-A*x
%
%or, for dense A, 
%
%   res = lssresidual(R,A,intval(speye(n)))
%
%Input   R     nxn real or complex point matrix
%        A     nxn real or complex point matrix
%Output  res   approximation of I-R*A
%
%Note that the fact that the last parameter is of type intval causes computation of an
%inclusion of the result rather than an approximation. Also note that only the last
%parameter is allowed to be of type intval, not one of the factors.
%
%This routine can be used in verifylss to improve accuracy of inclusion. Automatic use
%can be switched on and off, see intvalinit. Basically, the factors A and x are split
%into an upper and lower part and treated seperately. Fast way of splitting is taken from
%  T.J. Dekker: A floating-point technique for extending the available precision,
%    Numerische Mathematik 18:224-242, 1971.
%For timing and improvement of accuracy of inclusion see the table below.
%
%For randomly generated full matrices of dimension n and condition number 1e8, and 
%with b=a*(2*rand(n,1)-1) the following table lists the computing time t1 w/o and 
%t2 with improved residual in seconds on a 750 MHz Pentium III Laptop, as well as 
%the ratio t1/t2. Furthermore, the median of the radius of the inclusion is given 
%in r1 and r2 and the improvement r1/r2 in radius, respectively.
%
%     n     t1     t2      t2/t1     r1       r2      inclusion of b-A*x
%--------------------------------------------------
%   500   0.0018  0.0091    5.1   2.0e-16  1.2e-21
%  2000   0.013   0.117     8.8   3.3e-16  6.7e-21
%  5000   0.068   0.63      9.3   5.0e-16  5.3e-19
%
%Similar numbers (in seconds) for a sparse matrix of dimensions 1000, 5000 and 10000 with 
%some 20 nonzero elements per row are as follows:
%
%     n     t1     t2      t2/t1     r1       r2      inclusion of b-A*x
%--------------------------------------------------
%  5000   0.0021  0.016     7.6   9.6e-17  4.4e-22
% 20000   0.0062  0.057     9.2   9.1e-17  4.2e-22
% 50000   0.015   0.14      9.7   8.8e-17  4.0e-22
%
%Similar numbers for a random nxn matrix A, an approximate inverse R and timing for I-R*A:
%
%     n     t1     t2      t2/t1     r1       r2      inclusion of I-R*A
%--------------------------------------------------
%   500   0.054   0.17      3.1   4.2e-11  2.9e-16
%  2000   1.8     5.4       2.9   6.6e-12  8.5e-17
%  5000  26      79         3.0   1.6e-12  4.0e-15
%

% written  04/02/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 05/30/07     S.M. Rump  typo
% modified 08/27/12     S.M. Rump  complex part
% modified 10/16/12     S.M. Rump  comment
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isa(A,'intval')                  % A interval
    error('Input A must be no interval')
  end
  
  if isa(x,'intval')                  % x interval
    error('Input x must be no interval')
  end
  
  if isreal(A)                        % A real
    if isreal(x)                    % x real
      if b.complex                % b complex
        error('Input A and x is real, so complex b makes no sense')
      else                        % A,x,b real
        res.complex = 0;
        [res.inf,res.sup] = lssresid(A,x,b.inf,b.sup);    
        res.mid = [];
        res.rad = [];
      end
    else                            % A real, x complex
      if b.complex
        setround(-1)                % convert intval b to inf/sup
        binf = b.mid - b.rad;
        setround(1)
        bsup = b.mid + b.rad;
      else
        binf = b.inf;
        bsup = b.sup;
      end
      [resrealinf,resrealsup] = lssresid(A,real(x),real(binf),real(bsup));
      [resimaginf,resimagsup] = lssresid(A,imag(x),imag(binf),imag(bsup));
      res.complex = 1;
      res.inf = [];
      res.sup = [];
      setround(1)
      res.mid = complex( resrealinf + 0.5*(resrealsup-resrealinf) , resimaginf + 0.5*(resimagsup-resimaginf) );
      res.rad = mag(res.mid - complex( resrealinf , resimaginf ));
    end
  else                                % A complex, no interval
    if isreal(x)                    % x real
      if b.complex
        setround(-1)                % convert intval b to inf/sup
        binf = b.mid - b.rad;
        setround(1)
        bsup = b.mid + b.rad;
      else
        binf = b.inf;
        bsup = b.sup;
      end
      [resrealinf,resrealsup] = lssresid(real(A),x,real(binf),real(bsup));
      [resimaginf,resimagsup] = lssresid(imag(A),x,imag(binf),imag(bsup));
      res.complex = 1;
      res.inf = [];
      res.sup = [];
      setround(1)
      res.mid = complex( resrealinf + 0.5*(resrealsup-resrealinf) , resimaginf + 0.5*(resimagsup-resimaginf) );
      res.rad = mag(res.mid - complex( resrealinf , resimaginf ));
    else                            % A and x complex
      if b.complex
        setround(-1)                % convert intval b to inf/sup
        binf = b.mid - b.rad;
        setround(1)
        bsup = b.mid + b.rad;
      else
        binf = b.inf;
        bsup = b.sup;
      end
      [resrealinf,resrealsup] = lssresid([real(A) imag(A)],[real(x);-imag(x)],real(binf),real(bsup));
      [resimaginf,resimagsup] = lssresid([real(A) imag(A)],[imag(x);real(x)],imag(binf),imag(bsup));
      res.complex = 1;
      res.inf = [];
      res.sup = [];
      setround(1)
      res.mid = complex( resrealinf + 0.5*(resrealsup-resrealinf) , resimaginf + 0.5*(resimagsup-resimaginf) );
      res.rad = mag(res.mid - complex( resrealinf , resimaginf ));
    end
  end
  
  res = class(res,'intval');
  
  setround(rndold)
  

function [resinf,ressup] = lssresid(A,x,binf,bsup)
%Inclusion of the residual b-A*x, splits A and x and multiplies A*x in parts

setround(0)
factor = 68719476737;                   % heuristically optimal splitting 2^36+1

C = factor*A;
Abig = C - A;
A1 = C - Abig;                          % small (upper) part from A
A2 = A - A1;                            % A = A1+A2 exact splitting

x = -x;
y = factor*x;
xbig = y - x;
x1 = y - xbig;                          % small (upper) part from -x
x2 = x - x1;                            % x = x1+x2 exact splitting

setround(-1)
resinf = (A1*x1+binf)+(A1*x2+A2*x);
setround(1)
ressup = (A1*x1+bsup)+(A1*x2+A2*x);
