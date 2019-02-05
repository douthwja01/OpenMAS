function c = vertcat(varargin)
%VERTCAT      Implements  [a(1) ; a(2) ; ...]  for Taylor
%

% written  05/21/09     S.M. Rump
%

  a = taylor(varargin{1});
  c.size = a.size;
  index = reshape(1:prod(a.size),a.size)';
  c.t = a.t(:,index);             % transposed of first element a
  cols = c.size(2);

  for i=2:length(varargin)
    a = taylor(varargin{i});
    if cols~=a.size(2)
      error('dimension do not fit')
    end
    c.size(1) = c.size(1) + a.size(1);
    index = reshape(1:prod(a.size),a.size)';
    c.t = [ c.t a.t(:,index) ];   % horzcat for columnwise stored arrays
  end
  
  index = reshape(1:prod(c.size),fliplr(c.size))';
  c.t = c.t(:,index);             % transpose result

  c = class(c,'taylor');
