function r = power(a,b)
%POWER        Implements  a .^ b  for intervals
%
%For complex b, for complex a and non-integer b, and for components with negative real a 
%  and b containing non-integer, implementation by exp( b .* log(a) ).
%

% written  10/16/98     S.M. Rump
% modified 12/30/98     S.M. Rump  improved performance and even exponent
% modified 05/19/02     S.M. Rump  problem with vector b fixed (thanks to Arrigo)
% modified 08/06/02     S.M. Rump  complete redesign
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    NaN corrected
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 07/28/05     S.M. Rump  NaN corrected (thanks to John Pryce)
% modified 11/02/05     S.M. Rump  a^0 := 0 (thanks to Jörg Kubitz)
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 05/13/09     S.M. Rump  integer exponents, NaN^0
% modified 11/17/12     S.M. Rump  a^intval(b) correct (Thanks to T. Ogita)
%

  if isreal(b)
    if isa(b,'intval') & ( b.inf==b.sup )
      b = b.inf;
      a = intval(a);
    end
    if isa(b,'double') & isreal(b) & prod(size(b))==1 & b==round(b)
      if b==0
        r = typeadj( ones(size(a)) , typeof(a) );
        r(isnan(a)) = NaN;
      else                        % b is integer
        b_is_negative = b<0;
        % check b is even to ensure result is nonnegative
        if b==2*floor(b/2)
          b = b/2;
          b_is_even = 1;
        else
          b_is_even = 0;
        end
        b = abs(b) - 1;           % abs(b) is at least 1
        r = a;
        while b>0
          if mod(b,2)==1
            r = r.*a;
          end
          b = floor(b/2);
          if b~=0
            a = sqr(a);
          end
        end
        if b_is_even
          r = sqr(r);
        end
        if b_is_negative
          r = 1./r;
        end
      end
      return
    end
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end
  
  % treat complex a or b
  if ( ~isreal(a) ) | ( ~isreal(b) )
    r = exp( b .* log(intval(a)) );  
    setround(rndold)
    return
  end
  
  % adjust sizes
  if length(a)~=1
    if length(b)~=1
      if ~isequal(size(a),size(b))
        error('sizes of mantissa and exponent do not match')
      end
    else
      b = repmat(full(b),size(a));
    end
  else
    if length(b)~=1
      a = repmat(full(a),size(b));
    end
  end
  
  wng = warning;
  warning off
  
  if isa(b,'intval')                  % exponent b real interval
    b_int = ( b.inf==b.sup ) & ( b.inf==round(b.inf) );
    if isa(a,'intval')                % a is interval, b interval
      if any(b_int(:))
        [rinfb_int rsupb_int] = power1(a.inf(b_int),a.sup(b_int),b.inf(b_int));
      end
      aposbnonint = ( ~b_int ) & ( a.inf>=0 );
    else                              % a is non-interval, b interval
      if any(b_int(:))                % b integer
        [rinfb_int rsupb_int] = power1(a(b_int),a(b_int),b.inf(b_int));
      end
      aposbnonint = ( ~b_int ) & ( a>=0 );
    end
  else                                % exponent b real non-interval
    b_int = ( b==round(b) );          % integer exponents
    if isa(a,'intval')                % a is interval, b not interval
      if any(b_int(:))
        [rinfb_int rsupb_int] = power1(a.inf(b_int),a.sup(b_int),b(b_int));
      end
      aposbnonint = ( ~b_int ) & ( a.inf>=0 );
    else                              % a is non-interval, b not interval
      if any(b_int(:))                % b integer
        [rinfb_int rsupb_int] = power1(a(b_int),a(b_int),b(b_int));
      end
      aposbnonint = ( ~b_int ) & ( a>=0 );
    end
  end
  
  r = intval(zeros(size(a)));  
  
  if any(b_int(:))            
    %VVVV r(b_int) = infsup(rinfb_int,rsupb_int);
    s.type = '()'; s.subs = {b_int}; r = subsasgn(r,s,infsup(rinfb_int,rsupb_int));
    %AAAA Matlab V5.2 bug fix
  end
  
  if any(aposbnonint(:))              % real result
    %VVVV r(aposbnonint) = exp(b(aposbnonint).*log(a(aposbnonint)));
    s.type = '()'; s.subs = {aposbnonint}; r = subsasgn(r,s,exp(subsref(b,s).*log(subsref(a,s))));
    %AAAA Matlab V5.2 bug fix
  end
  
  index = ~b_int & ~aposbnonint;      % possible complex result
  if any(index(:))
    %VVVV r(index) = exp(b(index).*log(a(index)));
    s.type = '()'; s.subs = {index}; r = subsasgn(r,s,exp(subsref(b,s).*log(subsref(a,s))));
    %AAAA Matlab V5.2 bug fix
  end
  
  % take care of NaNs
  index = isnan(a) | isnan(b);
  if any(index(:))
    %VVVV r(index) = NaN;
    s.type = '()'; s.subs = {index}; r = subsasgn(r,s,NaN);
    %AAAA Matlab V5.2 bug fix
  end
  
  warning(wng)
    
  setround(rndold)


function [rinf,rsup] = power1(ainf,asup,b)
% interval a, integer b, size(a)=size(b)

rinf = ones(size(ainf));
rsup = rinf;
b_0 = ( b==0 );
if any(b_0(:))
  rinf(b_0) = 1;
  rsup(b_0) = 1;
end
b_even = ( mod(b,2)==0 ) & ~b_0;
if any(b_even(:))
  a_0 = b_even & ( ainf<0 ) & (asup>0 );
  if any(a_0)                         % 0 in a  &  b even
    rinf(a_0) = 0;
    [lb,ub] = powerint(-ainf(a_0),asup(a_0),abs(b(a_0)),1,1);
    rsup(a_0) = max(lb,ub);
  end    
  a_pos = b_even & (ainf>=0 );
  if any(a_pos(:))                       % a>=0  &  b even
    [rinf(a_pos) rsup(a_pos)] = powerint(ainf(a_pos),asup(a_pos),abs(b(a_pos)),-1,1);
  end
  a_neg = b_even & (asup<=0 );
  if any(a_neg(:))                       % a<=0  &  b even
    [rinf(a_neg) rsup(a_neg)] = powerint(-asup(a_neg),-ainf(a_neg),abs(b(a_neg)),-1,1);
  end
end
b_odd = ( mod(b,2)==1 );
if any(b_odd(:))
  a_pos = b_odd & ( ainf>=0 );
  if any(a_pos)
    [rinf(a_pos) rsup(a_pos)] = powerint(ainf(a_pos),asup(a_pos),abs(b(a_pos)),-1,1);
  end
  a_neg = b_odd & ( asup<=0 );
  if any(a_neg)             % careful with rounding
    [rinf(a_neg) rsup(a_neg)] = powerint(-ainf(a_neg),-asup(a_neg),abs(b(a_neg)),1,-1);
    rinf(a_neg) = -rinf(a_neg);
    rsup(a_neg) = -rsup(a_neg);
  end
  a_0 = b_odd & ~a_pos & ~a_neg;
  if any(a_0)               % careful with rounding
    [rinf(a_0) rsup(a_0)] = powerint(-ainf(a_0),asup(a_0),abs(b(a_0)),1,1);
    rinf(a_0) = -rinf(a_0);
  end
end
b_neg = ( b<0 );                      % treat negative b
if any(b_neg(:))
  indexinv = b_neg & ( rinf<=0 ) & ( rsup>=0 );
  setround(-1)
  dummy = 1./rsup(b_neg);
  setround(1)
  rsup(b_neg) = 1./rinf(b_neg);
  rinf(b_neg) = dummy;
  if any(indexinv(:))
    rinf(indexinv) = -inf;
    rsup(indexinv) = inf;
  end
end


function [r1,r2] = powerint(a1,a2,b,rnd1,rnd2)
% non-negative a1,a2, positive integer b, size(a_)=size(b), result ri=ai.^b with rounding rndi

r1 = a1;
r2 = a2;
b = b(:)-1;
while any(b>0)
  index = ( mod(b,2)==1 );
  if any(index)
    setround(rnd1)
    r1(index) = r1(index).*a1(index);
    setround(rnd2)
    r2(index) = r2(index).*a2(index);
  end
  b = floor(b/2);
  index = ( b~=0 );
  if any(index)
    setround(rnd1)
    a1(index) = a1(index).*a1(index);
    setround(rnd2)
    a2(index) = a2(index).*a2(index);
  end
end
