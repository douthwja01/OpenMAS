function r = isempty(c)
%ISEMPTY      Returns 1 if input is empty, i.e. [], in the sense of Matlab
%
%   r = isempty(a)
%
%There are no empty intervals in the mathematical sense in INTLAB. For details,
%  please see Readme.txt
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/12/05     S.M. Rump  functionality exchanged with @hessian\isempty_
% modified 04/22/09     S.M. Rump  isempty_ removed
%

  r = isempty(c.x) | isempty(c.dx)| isempty(c.hx);
