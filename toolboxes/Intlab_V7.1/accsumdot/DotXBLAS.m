function s = DotXBLAS(x,y)
%DOTXBLAS     Dot product 'as if' computed in 2-fold precision
%
%   res = DotXBLAS(x,y)
%
%On return, res approximates x'*y with accuracy as if computed 
%  in 2-fold precision.
%
%Implements algorithm BLAS_ddot_x from
%  X. Li, J. Demmel, D. Bailey, G. Henry, Y. Hida, J. Iskandar, 
%    W. Kahan, S. Kang, {S.}, A. Kapur, M. Martin, B. Thompson, {B.},
%    T. Tung, {T.}, D. Yoo: Design, Implementation and Testing of 
%    Extended and Mixed Precision BLAS, ACM Trans. Math. Software, 
%    2(28), p. 152-205, 2002.
%Requires 37n flops.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest, complex input
%

  if ( ~isreal(x) ) | ( ~isreal(y) )
    error('DotXBLAS for real input only')
  end
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  s = 0;
  t = 0;
  for i=1:length(x)
    [h,r] = TwoProduct(x(i),y(i));
    [s1,s2] = TwoSum(s,h);
    [t1,t2] = TwoSum(t,r);
    s2 = s2 + t1;
    [t1,s2] = FastTwoSum(s1,s2);
    t2 = t2 + s2;
    [s,t] = FastTwoSum(t1,t2);
  end
  
  if rndold
    setround(rndold)
  end
  