function u = vertcat(varargin)
%VERTCAT      Implements  [a(1) ; a(2) ; ...]  for long numbers
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  u = long(varargin{1});

  for i=2:length(varargin)
    a = long(varargin{i});

    % if necessary, adapt length of mantissa
    preca = size(a.mantissa,2);
    precu = size(u.mantissa,2);
    if preca<precu
      a.mantissa(end,precu) = 0;
    elseif preca>precu
      u.mantissa(end,preca) = 0;
    end

    u.sign = [ u.sign ; a.sign ];
    u.mantissa = [ u.mantissa ; a.mantissa ];
    u.exponent = [ u.exponent ; a.exponent ];
    u.error.mant = [ u.error.mant ; a.error.mant ];
    u.error.exp = [ u.error.exp ; a.error.exp ];
  end
