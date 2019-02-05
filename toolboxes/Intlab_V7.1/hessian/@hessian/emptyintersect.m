function [empty,a] = emptyintersect(a,b)
%EMPTYINTERSECT    Compute intersection and check for empty components
%
%   [empty,c] = emptyintersect(a,b)
%
%Result c is that of intersect(a,b), and 
%  empty(i) = 1     intersection of a(i) and b(i) is empty
%             0     intersection of a(i) and b(i) is not empty
%             NaN   at least one of a(i) and b(i) is NaN
%
%Intersection is taken for a.x and b.x part, derivatives remain unchanged.
%Input a and b must be both real or both complex
%

% written  04/23/12     S.M. Rump  (Thanks to Kolumbán Sándor for pointing
%                                     to that missing function)
%

  a = hessian(a);
  b = hessian(b);

  if isa(a.x,'intval') | isa(b.x,'intval')
    a = intval(a);
  end

  [empty,a.x] = emptyintersect(a.x,b.x);
  a.dx(:,find(empty)) = NaN;
  a.hx(:,find(empty)) = NaN;
  