function c = subsref(a,s)
%SUBSREF      Implements subscripted references for long numbers
%
%   example   c = a(3:5)
%
%Only one index because long array is always column vector
%

%Access to subfields sign, exponent, mantissa, error for internal use only
%

% written  12/30/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  INTLAB_LONG_ERROR = getappdata(0,'INTLAB_LONG_ERROR');

  while 1
    if strcmp(s(1).type,'()')     % index reference a(i)
      if length(s(1).subs)>1
        error('multiple index in long subsref')
      end
      if isa(a,'double')              % for a.sign(i) etc.
        c = a(s(1).subs{:});
      else
        c.sign = a.sign(s(1).subs{:});
        c.exponent = a.exponent(s(1).subs{:});
        c.mantissa = a.mantissa(s(1).subs{:},:);
        if INTLAB_LONG_ERROR
          c.error.mant = a.error.mant(s(1).subs{:});
          c.error.exp = a.error.exp(s(1).subs{:});
        else
          c.error.mant = 0;
          c.error.exp = 0;
        end
        c = class(c,'long');
      end
    elseif strcmp(s(1).type,'.')     % subfield access
      if     strcmp(s(1).subs,'sign'), c = a.sign;
      elseif strcmp(s(1).subs,'exponent'), c = a.exponent;
      elseif strcmp(s(1).subs,'mantissa'), c = a.mantissa;
      elseif strcmp(s(1).subs,'error') & INTLAB_LONG_ERROR, c = a.error;
      else
        error('invalid subscript reference for intval')
      end
    else
      error('invalid index reference for intval')
    end
    if length(s)==1  
      if rndold
        setround(rndold)
      end
      return
    end
    s = s(2:end);
    a = c;
  end
