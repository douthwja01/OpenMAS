function res = isign(x)
%ISIGN        Interval sign, internal function for slope abs
%

%sign of interval array:
%  +1      for x>=0
%  -1      for x<=0
%  [-1,1]  for 0 in x

% written  12/06/98     S.M. Rump
% modified 12/27/02     S.M. Rump  resinf,ressup changed to type double: fix of Matlab 6.5 discrepancy to 6.0
%

  xinf = inf(x);
  xsup = sup(x);
  resinf = double( xinf>=0 );
  resinf(xsup<=0) = -1;
  ressup = resinf;
  index = ( xinf<0 ) & ( xsup>0 );
  if any(index(:))
    resinf(index) = -1;
    ressup(index) = 1;
  end
  res = infsup(resinf,ressup);
