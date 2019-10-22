function r = isempty(c)
%ISEMPTY      Returns 1 if input is empty, i.e. [], in the sense of Matlab
%
%   r = isempty(a)
%
%There are no empty intervals in the mathematical sense in INTLAB. For details,
%  please see Readme.txt
%

% written  05/21/09     S.M. Rump
%

  r = isempty(c.t);
