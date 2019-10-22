function rnd = getround
%GETROUND     Get rounding mode
%
%   rnd = getround
%
%On return,
%
%   rnd = -1   for rounding downwards
%   rnd =  0   for rounding to nearest
%   rnd =  1   for rounding upwards
%   rnd =  2   for rounding towards zero (chop)
%

% written  11/23/98     S.M. Rump
% modified 12/04/05     S.M. Rump  improved performance
% modified 12/15/07     T. Ogita   modified for Intel-based Mac
%
 
e = [];
e = 1e-30;
x = 1 + e;
y = 1 - e;
if x == y                      % fast check for rounding to nearest
    rnd = 0;
else
    e = [];
    e = 1e-30;
    z = (-1) + e;
    if ( x==1 ) & ( z==-1 )    % round downwards
        rnd = -1;
    elseif ( y==1 )            % round upwards
        rnd = +1;
    else
        rnd = 2;
    end
end
