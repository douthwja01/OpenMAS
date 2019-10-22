function res = relacc(x)
%RELACC       Relative accuracy of an interval in decimal digits
%
%   res = relacc(x)   for interval input x
%
%Relative accuracy as defined in
%  S.M. Rump: The relative accuracy of an interval, to appear
%In numerical analysis we have the well-accepted rule of thumb that a computed approximation 
%  of a linear system with condition number 10^k has roughly 16-k correct digits.
%The purpose of relacc is to check this rule of thumb for interval inclusions. That means 
%  a good verification algorithm produces an inclusion x with relacc(x) roughly of size 16-k.
%
%Note that r(i)=relacc(x(i)) for r=relacc(x) is only true for normed vector x, more precisely,
%  for max(rad(x)+mig(x))==1. For real x this means max(abs(inf(x))+abs(sup(x)))==1.
%
%Careful: relacc produces full output for sparse x; use relacc(nonzeros(x)) for sparse input
%

% written  10/20/12     S.M. Rump

  if x.complex                  % complex interval input
        
    if isequal(x.rad,0)         % point interval
      
      res = inf(size(x.mid));   % infinitely accurate
      res(isnan(x.mid)) = NaN;  % treat NaN midpoint separately
      
    else
      
      N = x.rad + mig(x);       % works for vectors as well
      N = max( max(N(:)) , 1 );
      res = -log10( x.rad ./ N );
      
    end
        
  else                          % real interval input
    
    N = abs(x.inf) + abs(x.sup);
    N = 2*max( max(N(:)) , 1 );
    res = -log10( ( x.sup - x.inf ) ./ N );
    
end
    