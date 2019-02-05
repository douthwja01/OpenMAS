function a = band(a,p,q)
%BAND         Extract band from matrix a, lower bandwidth p, upper bandwidth q
%   if parameter q is omitted, q:=p
%
%   c = band(a,p,q)
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin<3
    q = p;
  end

  index = ( tril(triu(ones(size(a.x)),-p),q) == 0 );
  a.x(index) = 0;
  a.dx(:,index) = 0;
  a.hx(:,index) = 0;
