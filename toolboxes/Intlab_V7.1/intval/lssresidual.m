function res = lssresidual(A,x,b)
%LSSRESIDUAL    Improved approximation of residual b-A*x or I-R*A (poor men's residual)
%
%Typical calls are
%
%   res = lssresidual(A,x,b)
%
%Input   A     nxn real or complex point matrix (dense or sparse)
%        x     some approximation to A\b
%        b     real or complex point vector or nxk matrix
%Output  res   approximation of b-A*x
%
%or, for dense or sparse A, 
%
%   res = lssresidual(R,A,speye(n))
%
%Input   R     nxn real or complex point matrix
%        A     nxn real or complex point matrix
%Output  res   approximation of I-R*A
%
%This routine can be used in verifylss to improve accuracy of inclusion. Automatic use
%can be switched on and off, see intvalinit. Basically, the factors A and x are split
%into an upper and lower part and treated separately. Fast way of splitting is taken from
%  T.J. Dekker: A floating-point technique for extending the available precision,
%    Numerische Mathematik 18:224-242, 1971.
%For timing and improvement of accuracy of inclusion see the table below.
%
%For randomly generated full matrices of dimension n with b=a*randn(n,1) the following
%table lists the computing time t1 w/o and t2 with improved residual in seconds on a 
%2.8 GHz Laptop, as well as the ratio t2/t1.
%
%     n     t1      t2    t2/t1     approximation of b-A*x
%-------------------------------
%   500   0.0001  0.0069   59.1
%  2000   0.0044  0.097    21.9 
%  5000   0.028   0.54     19.3  
%
%Similar numbers (in seconds) for a sparse matrix of dimensions 1000, 5000 and 10000 with 
%some 20 nonzero elements per row are as follows:
%
%     n     t1      t2    t2/t1     approximation of b-A*x
%-------------------------------
%  5000   0.0003  0.013    42.7
% 20000   0.0018  0.051    27.9
% 50000   0.0063  0.13     20.0  
%
%Similar numbers for a random nxn matrix A, an approximate inverse R and timing for I-R*A:
%
%     n     t1      t2    t2/t1     approximation of I-R*A
%-------------------------------
%   500   0.017   0.086     4.9
%  2000   0.92    2.8       3.1   
%  5000  13      40         3.1
%

% written  04/02/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 12/12/05     S.M. Rump  splitting corrected
% modified 05/30/07     S.M. Rump  typo
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/16/12     S.M. Rump  comment
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if isreal(A)                          % A real
    if isreal(x)                        % x real
      if isreal(b)                      % b real
        res = resid(A,x,b);    
      else
        error('Input A and x is real, so complex b makes no sense')
      end
    else                                % A real, x complex
      res = complex( resid(A,real(x),real(b)) , resid(A,imag(x),imag(b)) );
    end
  else                                  % A complex
    if isreal(x)                        % x real
      res = complex( resid(real(A),x,real(b)) , resid(imag(A),x,imag(b)) );
    else                                % A and x complex
      res = complex( resid([real(A) imag(A)],[real(x);-imag(x)],real(b)) , resid([real(A) imag(A)],[imag(x);real(x)],imag(b)) );
    end
  end
  
  if rndold
    setround(rndold)
  end


function res = resid(A,x,b)
%Approximation of the residual b-A*x, splits A and x and multiplies A*x in parts
%

  setround(0)
  factor = 68719476737;                   % heuristically optimal splitting 2^36+1

  C = factor*A;
  Abig = C - A;
  A1 = C - Abig;                          % small (upper) part of A
  A2 = A - A1;                            % A = A1+A2 exact splitting

  x = -x;
  y = factor*x;
  xbig = y - x;
  x1 = y - xbig;                          % small (upper) part of -x
  x2 = x - x1;                            % -x = x1+x2 exact splitting

  res = (A1*x1+b)+(A1*x2+A2*x);
