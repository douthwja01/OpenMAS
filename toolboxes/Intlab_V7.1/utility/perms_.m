function P = perms_(x)
%PERMS_       Functionality like MATLAB\PERMS, but fast
%
%   P = perms_(x)
%
% permuations of x in array of size n! x n for length(x)=n
%
%Sample timing comparing perms_(1:n) with built-in Matlab routine perms(1:n)
%in seconds on 120 MHz Pentium Laptop:
%
%
%      n     t(perms)  t(perms_)
% --------------------------------
%   2.0000    0.0052    0.0044
%   3.0000    0.0198    0.0091
%   4.0000    0.0875    0.0154
%   5.0000    0.4390    0.0253
%   6.0000    2.6650    0.0495
%   7.0000   18.8400    0.1870
%   8.0000  167.4100    1.3400
%
%
%Matlab description:
%   PERMS(1:N), or PERMS(V) where V is a vector of length N, creates a
%   matrix with N! rows and N columns containing all possible
%   permutations of the N elements.
%
%   This function is only practical for situations where N is less
%   than about 15.
%
%   See also NCHOOSEK, RANDPERM, PERMUTE.


% written  03/19/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/30/07     S.M. Rump  typo
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  n=length(x);
  if n==0
    P = [];
  else
    P = zeros(prod(1:n),n);
    P(1,1)=1;
    fac = 1;
    for k=2:n
      Pold = P(1:fac,1:k-1);
      w = 1:k;
      v = ((k-1)*fac+1):(k*fac);
      for i=k:-1:1
        P(v,w) = [ i*ones(fac,1) Pold+(Pold>=i) ];
        v = v-fac;
      end
      fac = k*fac;
    end
    P = x(P);
  end
  
  if rndold
    setround(rndold)
  end
