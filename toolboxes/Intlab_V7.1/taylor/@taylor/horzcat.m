function c = horzcat(varargin)
%HORZCAT      Implements  [a(1) a(2) ...]  for Taylor
%

% written  05/21/09     S.M. Rump
%

  a = taylor(varargin{1});
  c.size = a.size;
  c.t = a.t;
  rows = c.size(1);

  for i=2:length(varargin)
    a = taylor(varargin{i});
    if rows~=a.size(1)
      error('dimension do not fit')
    end
    c.size(2) = c.size(2) + a.size(2);
    c.t = [ c.t a.t ];     % uses that arrays are stored columnwise
  end

  c = class(c,'taylor');
