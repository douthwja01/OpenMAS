function E = exp2index(e,emax)
%EXP2INDEX    internal routine for polynom
%
%For an m x n array e of exponents, 0<=e(:,i)<emax(i) for row n-vector emax,
%  output E is a column m-vector of linearized indices, 0 <= E(:) < prod(emax)
%
  m = size(e,1);
  n = length(emax);
  E = sum( e .* [ ones(m,1) cumprod(ones(m,1)*emax(1:n-1),2) ] , 2 );
