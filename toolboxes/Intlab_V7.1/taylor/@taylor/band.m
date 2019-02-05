function a = band(a,p,q)
%BAND         Extract band from matrix a, lower bandwidth p, upper bandwidth q
%   if parameter q is omitted, q:=p
%
%   c = band(a,p,q)
%

% written  05/21/09     S.M. Rump
%

  if nargin<3
    q = p;
  end

  index = reshape(1:prod(a.size),a.size);
  index = find(tril(triu(index,-p),q)==0);
  a.t(:,index) = 0;
