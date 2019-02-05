function r = isempty(a)
%ISEMPTY      Returns 1 if input is empty, i.e. [], in the sense of Matlab
%
%   r = isempty(a)
%
%There are no empty intervals in the mathematical sense in INTLAB. For details,
%  please see Readme.txt
%

% written  10/15/99     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 10/12/05     S.M. Rump  functionality exchanged with @intval\isempty_
% modified 04/22/09     S.M. Rump  isempty_ removed
%

  if a.complex
    r = isempty(a.mid) | isempty(a.rad);
  else
    r = isempty(a.inf) | isempty(a.sup);
  end
