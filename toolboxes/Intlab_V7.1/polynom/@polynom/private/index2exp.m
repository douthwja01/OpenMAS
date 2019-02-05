function e = index2exp(E,emax)
%INDEX2EXP    internal routine for polynom
%
%Reverses exp2index, result index to E w.r.t. emax
%

  m = length(E);
  n = length(emax);
  e = zeros(m,n);
  b = [ 1 cumprod(emax(1:n-1)) ];
  for i=n:-1:1
    e(:,i) = floor(E/b(i));
    E = rem(E,b(i));
  end
