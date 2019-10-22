function c = horzcat(varargin)
%HORZCAT      Implements  [a(1) a(2) ...]  for gradients
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  a = gradient(varargin{1});
  c.x = a.x;
  c.dx = a.dx;

  for i=2:length(varargin)
    a = gradient(varargin{i});
    c.x = [ c.x , a.x ];
    c.dx = [ c.dx ; a.dx ];     % uses that arrays are stored columnwise
  end

  c = class(c,'gradient');
