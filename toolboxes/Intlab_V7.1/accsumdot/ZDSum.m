function y = ZDsum(v)
%ZDSUM        Computes approximation of sum(p) using Zielke/Drygalla's approach
%
%   y = ZDsum(p)
%
%On return, res is an approximation of sum(p) of the column vector p
%  provided no over- or underflow occurred. 
%
%Taken from Algorithm AccDot on page 30 in
%  G. Zielke, V. Drygalla: Genaue Lösung linearer Gleichungssysteme, 
%    GAMM Mitteilungen 26(1-2), p. 8-107, 2003.
%Requires (6m+3)n+7k flops for m executions of while-loop and
%  length(s)=k for intermediate vector s.
%
%Reference implementation! Slow due to interpretation!
%

% written  03/03/07     S.M. Rump
% modified 05/09/09     S.M. Rump  rounding to nearest
%
  
  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  n = length(v);
  if isa(v,'double')
    C1 = 54;
    C2 = -1023;
  else
    C1 = 25;
    C2 = -127;
  end
  ma = max(abs(v));
  [mmax,emax] = log2(ma);
  q = 2^emax;
  k = floor( C1 - log(n)/log(2) );
  p = 2^k;
  i = 0;
  while any(v ~= zeros(n,1)) & ( q/p > 2^C2 )
    i = i+1;
    q = q/p;
    g = fix(v/q);
    s(i) = sum(g);
    v = v - g*q;
  end
  i = i+1;
  s(i) = sum(v)*p/q;
  ue = 0;
  for j=i:(-1):1
    t = s(j) + ue;
    ue = floor(t/p);
    s(j) = t - ue*p;
  end
  y = ue;
  for j=1:i
    y = s(j) + y*p;
  end
  y = y*q/p;
  
  if rndold
    setround(rndold)
  end
