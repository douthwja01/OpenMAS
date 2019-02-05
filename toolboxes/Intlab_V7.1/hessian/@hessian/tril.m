function a = tril(a,k)
%TRIL         Implements  tril(a,k)  for hessians
%
%   c = tril(a,k)
%
% functionality as Matlab function tril for matrices
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin==1
    k = 0;
  end

  a.x = tril(a.x,k);
  index = ( tril( ones(size(a.x)) , k ) == 0 );
  a.dx(:,index) = 0;
  a.hx(:,index) = 0;
