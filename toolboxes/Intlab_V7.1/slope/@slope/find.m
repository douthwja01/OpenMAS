function [ia,ja,sa] = find(a)
%FIND         Implements  find(a)  for slopes
%
%   index = find(a)
%   [ia,ja] = find(a)
%   [ia,ja,sa] = find(a)
%
% functionality as Matlab function find
%

% written  12/06/05     S.M. Rump
%

  wng = warning;
  warning off

  if nargout<=1
    ia = find( reshape( any(a.r,2) | any(a.s,2) , a.size ) );
  elseif nargout==2
    [ia,ja] = find( reshape( any(a.r,2) | any(a.s,2) , a.size ) );
  else
    ia = find( reshape( any(a.r,2) | any(a.s,2) , a.size ) );
    %VVVV  sa = a(ia);
    s.type = '()'; s.subs = {ia}; sa = subsref(a,s);
    %AAAA  Matlab bug fix
    [ia,ja] = find( reshape( any(a.r,2) | any(a.s,2) , a.size ) );
  end

  warning(wng)
  