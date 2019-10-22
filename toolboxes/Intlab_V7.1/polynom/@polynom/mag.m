function p = mag(p);
%MAG          Polynomial with absolute values of coefficients
%
%   r = mag(p);
%
%r_i = mag(p_i), i.e. result is real polynomial, also for interval input
%Replaces abss. Thanks to Arnold for pointing to better name.
%

% written  01/04/05     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  p.c = mag(p.c);
