function r = isempty(a)
%ISEMPTY      Returns 1 if input is empty, i.e. [], in the sense of Matlab
%
%   r = isempty(a)
%
%There are no empty intervals in the mathematical sense in INTLAB. For details,
%  please see Readme.txt
%

% written  09/29/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/12/05     S.M. Rump  functionality exchanged with @long\isempty_
% modified 04/22/09     S.M. Rump  isempty_ removed
%

  r = isempty(a.sign) | isempty(a.mantissa) | isempty(a.exponent);
