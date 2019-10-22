function a = band(a,p,q)
%BAND         Extract band from matrix a, lower bandwidth p, upper bandwidth q
%   if parameter q is omitted, q:=p
%
%   c = band(a,p,q)
%

% written  09/28/01     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    improved performance
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin<3
    q = p;
  end

  a = tril(triu(a,-p),q);
