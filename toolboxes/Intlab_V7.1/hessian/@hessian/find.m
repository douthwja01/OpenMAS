function [ia,ja,sa] = find(a)
%FIND         Implements  find(a)  for hessians
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
    ia = union( find(a.x) , find( any(a.dx,1) | any(a.hx,1) ) );
  elseif nargout==2
    [ia,ja] = find( logical(a.x) | reshape( any(a.dx,1) | any(a.hx,1) , size(a.x) ) );
  else
    ia = union( find(a.x) , find( any(a.dx,1) | any(a.hx,1) ) );
    %VVVV  sa = a(ia);
    s.type = '()'; s.subs = {ia}; sa = subsref(a,s);
    %AAAA  Matlab bug fix
    [ia,ja] = find( logical(a.x) | reshape( any(a.dx,1) | any(a.hx,1) , size(a.x) ) );
  end

  warning(wng)
  