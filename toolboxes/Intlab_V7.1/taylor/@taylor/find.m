function [ia,ja,sa] = find(a)
%FIND         Implements  find(a)  for Taylor
%
%   index = find(a)
%   [ia,ja] = find(a)
%   [ia,ja,sa] = find(a)
%
% functionality as Matlab function find
%

% written  05/21/09     S.M. Rump
%

  if nargout<=1
    ia = find(any(a.t));
  elseif nargout==2
    index = reshape(1:prod(a.size),a.size);
    index(~any(a.t)) = 0;
    [ia,ja] = find(index);
  else
    ia = find(any(a.t));
    %VVVV  sa = a(ia);
    s.type = '()'; s.subs = {ia}; sa = subsref(a,s);
    %AAAA  Matlab bug fix
    index = reshape(1:prod(a.size),a.size);
    index(~any(a.t)) = 0;
    [ia,ja] = find(index);
  end
