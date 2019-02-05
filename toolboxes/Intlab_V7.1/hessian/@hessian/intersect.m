function a = intersect(a,b)
%INTERSECT    Intersection of Hessians; empty components set to NaN
%
%   c = intersect(a,b)
%
%Result is an ***outer*** inclusion of the a.x and b.x components, the 
%  derivates remain unchanged. For details, see intval/intersect.
%
%===> To check empty intersection please use the INTLAB function
%===>   "EmptyIntersect" 
%
%Input a and b must be both real or both complex
%

% written  04/23/12     S.M. Rump  (Thanks to Kolumbán Sándor for pointing
%                                     to that missing function)
%

  [empty,a] = emptyintersect(a,b);
