function Pi = longpi
%Calculate Pi to current longprecision
%
%Typical call:
%   longprecision(200), p=longpi, display(p,0)
%

% written  11/30/98     S.M. Rump
% modified 01/10/04     S.M. Rump  faster approach for large n added
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
%

  if 0  % approximation by sine formula
    INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
    L = longprecision;
    k = 10;
    a = long(pi);
    while k<L/2
      k = 2*k;
      if k<L
        longinit('WithoutErrorTerm',0);
      else
        longinit('WithErrorTerm',0);
      end
      % longprecision(k);
      % a = a + sin(a);
      % y = sin(a)
      t = a;
      aa = a*a;
      apos = longshift(a,1);
      aneg = 0;
      i = 1;
      while 1
        t = t*aa/(2*i*(2*i+1));
        if mod(i,2)==0
          apos = apos + t;
        else
          aneg = aneg + t;
        end
        if t.exponent < -INTLAB_LONG_PRECISION-1
          a = apos - aneg;
          break
        end
        i = i+1;
      end
    end
    Pi = a;
    return
  end

  if longprecision>1000
    % A standard approach
    Pi = 48*atan1(18) + 32*atan1(57) - 20*atan1(239) ;
  else
    % A faster formula due to Fabrice Bellard:
    INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');
    F = inline('      ( (((((long(16400000*n+60500000))*n+90764000)*n+70652400)*n+29980024)*n+6543234)*n+570042 )/      ( long(5*(4*n+1)*(4*n+3))*long((10*n+1)*(10*n+3))*long((2*n+1)*(10*n+7))*(10*n+9) ) ');
    s = 0;
    n = 0;
    t = long(1);
    while ( t.exponent > -INTLAB_LONG_PRECISION-1 )
      t = longshift(F(n),-10*n);
      if odd(n)
        s = s-t;
      else
        s = s+t;
      end
      n = n+1;
    end
    Pi = longshift(s,-6);
  end


  function s = atan1(x)
    % y = arctan(1/x)
    INTLAB_LONG_PRECISION = getappdata(0,'INTLAB_LONG_PRECISION');

  t = long(1)/x;
  spos = t;
  sneg = 0;
  i = 1;
  while 1
    t = t/( x*x );
    if mod(i,2)==0
      spos = spos + t/(2*i+1);
    else
      sneg = sneg + t/(2*i+1);
    end
    if t.exponent < -INTLAB_LONG_PRECISION-1
      s = spos - sneg;
      break
    end
    i = i+1;
  end
  