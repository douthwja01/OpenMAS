function c = subsasgn(c,s,b)
%SUBSASGN     Implements subscripted assignments for long numbers
%
%   example   c(2:3) = b
%

% written  12/30/98     S.M. Rump
% modified 11/11/99     S.M. Rump  precb<=precc
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  if length(s)>1
    error('multiple indexing for long assignment not allowed')
  end

  if strcmp(s.type,'()')     % assignment c(i) = b
    if isempty(b)
      c.sign(s(1).subs{:}) = [];
      c.exponent(s(1).subs{:}) = [];
      c.mantissa(s(1).subs{:},:) = [];
      if INTLAB_LONG_ERROR
        c.error.mant(s(1).subs{:}) = [];
        c.error.exp(s(1).subs{:}) = [];
      end
      return
    else
      if isa(b,'double')
        b = long(b);
      end
    end
    if isempty(c)
      c.sign = [];
      c.exponent = [];
      c.mantissa = [];
      c.error.mant = [];
      c.error.exp = [];
      c = class(c,'long');
    end
    c.sign(s(1).subs{:},1) = b.sign;
    c.exponent(s(1).subs{:},1) = b.exponent;
    precb = size(b.mantissa,2);
    precc = size(c.mantissa,2);
    if precb<=precc
      c.mantissa(s(1).subs{:},precb+1:end) = 0;
      c.mantissa(s(1).subs{:},1:precb) = b.mantissa;
    else
      c.mantissa(end,precb) = 0;
      c.mantissa(s(1).subs{:},:) = b.mantissa;
    end
    if INTLAB_LONG_ERROR
      c.error.mant(s(1).subs{:},1) = b.error.mant;
      c.error.exp(s(1).subs{:},1) = b.error.exp;
    end
  else
    error('invalid call of subsasgn')
  end
