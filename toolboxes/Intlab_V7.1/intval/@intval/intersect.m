function c = intersect(a,b)
%INTERSECT    Intersection of intervals; empty components set to NaN
%
%   c = intersect(a,b)
%
%Result is an ***outer*** inclusion:
%   c includes the true intersection of a and b
%   components of c are set to NaN, if the corresponding components of
%     a and b are NaN, or if the intersection is empty
%   Thus a component NaN of the result c does NOT mean an empty intersection,
%     but rather that no information is available.
%
%===> To check empty intersection please use the INTLAB function
%===>   "EmptyIntersect" 
%
%Input a and b must be both real or both complex
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/06/05     S.M. Rump  correction of result (thanks to Andreas Rauh, Ulm, for hint)
% modified 06/11/06     S.M. Rump  correction of result (thanks to Andreas Rauh, Ulm, for hint)
% modified 09/10/07     S.M. Rump  redesign
% modified 04/22/09     S.M. Rump  comment changed, based on emptyintersect
%

  [empty,c] = emptyintersect(a,b);
