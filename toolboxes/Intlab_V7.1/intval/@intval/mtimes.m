function c = mtimes(a,b,dummy)
%MTIMES       Implements  a * b  for intervals
%

% written  10/16/98     S.M. Rump
% modified 12/30/98     S.M. Rump  improved performance
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    major revision
%                                    remove check for 'double', sparse input
%                                    extra argument for non-interval output (internal use only)
%                                    improved performance for outer products
%                                    take care of Matlab sparse Inf/NaN bug
% modified 08/24/04     S.M. Rump  IR x IR for sharp multiplication 
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
% modified 05/23/06     S.M. Rump  sparse Inf/NaN bug corrected in Version 7.2+
% modified 12/03/06     S.M. Rump  Sparse Bug global flag (thanks to Arnold)
% modified 12/05/06     S.M. Rump  Fast sharp multiplication if one factor contains no zero-IVs
% modified 05/09/09     S.M. Rump  thin intervals with equal NaNs
% modified 08/20/12     S.M. Rump  performance improvement, in particular for SharpIVMult 
%                                    (thanks to Adam Kleiner for pointing to this), and 
%                                    access to sparse submatrices
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_INTVAL_IVMULT = getappdata(0,'INTLAB_INTVAL_IVMULT');

  [m n] = size(a); 
  [n1 n2] = size(b);

  if ( m*n==1 ) | ( n1*n2==1 )
    if nargin==2
      c = a .* b;
    else
      c = times(a,b,0);
    end
    return
  end

  if n~=n1
    error('matrices not compatible for multiplication')
  end

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
  end

  if ~isa(a,'intval')                     % a is double
    if ~isreal(a) | b.complex             % complex case
      if ~b.complex
        setround(1)
        b.mid = b.inf + 0.5*(b.sup-b.inf);
        b.rad = b.mid - b.inf;
      end
      c.complex = 1;                      % R * IC  or  C * IC
      c.inf = [];
      c.sup = [];
      if isreal(a) | ~b.complex           % one real factor
        setround(-1)                      % R * IC
        c1 = a * b.mid;
        setround(1)
        c2 = a * b.mid;
      else                                % C * IC
        setround(-1)
        c1 = real(a) * real(b.mid) + (-imag(a)) * imag(b.mid) + ...
             ( real(a) * imag(b.mid) + imag(a) * real(b.mid) ) * j;
        setround(1)
        c2 = real(a) * real(b.mid) + (-imag(a)) * imag(b.mid) + ...
             ( real(a) * imag(b.mid) + imag(a) * real(b.mid) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % R * IC  or  C * IC
      if isequal(b.rad,0)
        c.rad = abs(c.mid-c1);
        if ~any(find(c.rad))              % take care of huge arrays
          c.rad = 0;
        end
      else
        c.rad = abs(c.mid-c1) + abs(a) * b.rad;
      end
    else                                  % real case  R * IR
      c.complex = 0;
      bthin = isequalwithequalnans(b.inf,b.sup);
      if bthin                            % R * R with directed rounding
        setround(-1)
        c.inf = a*b.inf;
        setround(1)
        c.sup = a*b.inf;
      else
        setround(1)                       % R * IR
        bmid = b.inf + 0.5*(b.sup-b.inf);
        brad = bmid - b.inf;
        crad = abs(a) * brad;
        c.inf = 0;                        % preserve order of definition of intval c
        c.sup = a * bmid + crad;
        setround(-1)
        c.inf = a * bmid - crad;
      end
      c.mid = [];
      c.rad = [];
    end
  elseif ~isa(b,'intval')                 % b is double
    if a.complex | ~isreal(b)             % complex case
      if ~a.complex
        setround(1)
        a.mid = a.inf + 0.5*(a.sup-a.inf);
        a.rad = a.mid - a.inf;
      end
      c.complex = 1;                      % IC * R  or IC * C
      c.inf = [];
      c.sup = [];
      if ~a.complex | isreal(b)           % one real factor
        setround(-1)                      % IC * R
        c1 = a.mid * b;
        setround(1)
        c2 = a.mid * b;
      else                                % IC * C
        setround(-1)
        c1 = real(a.mid) * real(b) + (-imag(a.mid)) * imag(b) + ...
             ( real(a.mid) * imag(b) + imag(a.mid) * real(b) ) * j;
        setround(1)
        c2 = real(a.mid) * real(b) + (-imag(a.mid)) * imag(b) + ...
             ( real(a.mid) * imag(b) + imag(a.mid) * real(b) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % IC * R  or  IC * C
      if isequal(a.rad,0)
        c.rad = abs(c.mid-c1);
        if ~any(find(c.rad))              % take care of huge arrays
          c.rad = 0;
        end
      else
        c.rad = abs(c.mid-c1) + a.rad * abs(b);
      end
    else                                  % real case  IR * R
      c.complex = 0;
      athin = isequalwithequalnans(a.inf,a.sup);
      if athin                            % R * R with directed rounding
        setround(-1)
        c.inf = a.inf*b;
        setround(1)
        c.sup = a.inf*b;
      else
        setround(1)                       % IR * R
        amid = a.inf + 0.5*(a.sup-a.inf);
        arad = amid - a.inf;
        crad = arad * abs(b);
        c.inf = 0;                        % preserve order of definition of intval c
        c.sup = amid * b + crad;
        setround(-1)
        c.inf = amid * b - crad;
      end
      c.mid = [];
      c.rad = [];
    end
  else                                    % both a and b interval
    if a.complex | b.complex              % complex case
      if ~a.complex
        setround(1)
        a.mid = a.inf + 0.5*(a.sup-a.inf);
        a.rad = a.mid - a.inf;
      end
      if ~b.complex
        setround(1)
        b.mid = b.inf + 0.5*(b.sup-b.inf);
        b.rad = b.mid - b.inf;
      end
      c.complex = 1;
      c.inf = [];
      c.sup = [];
      if ~a.complex | ~b.complex          % one real factor
        setround(-1)
        c1 = a.mid * b.mid;
        setround(1)
        c2 = a.mid * b.mid;
      else                                % IC * IC
        setround(-1)
        c1 = real(a.mid) * real(b.mid) + (-imag(a.mid)) * imag(b.mid) + ...
             ( real(a.mid) * imag(b.mid) + imag(a.mid) * real(b.mid) ) * j;
        setround(1)
        c2 = real(a.mid) * real(b.mid) + (-imag(a.mid)) * imag(b.mid) + ...
             ( real(a.mid) * imag(b.mid) + imag(a.mid) * real(b.mid) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % IC * IC
      if isequal(a.rad,0)
        if isequal(b.rad,0)
          c.rad = abs(c.mid-c1);
          if ~any(find(c.rad))            % take care of huge arrays
            c.rad = 0;
          end
        else
          c.rad = abs(c.mid-c1) + abs(a.mid) * b.rad;
        end
      elseif isequal(b.rad,0)
        c.rad = abs(c.mid-c1) + a.rad * abs(b.mid);
      else
        c.rad = abs(c.mid-c1) + ...
                  a.rad * ( abs(b.mid) + b.rad ) + abs(a.mid) * b.rad;
      end
    else                                  % real case,  IR * IR
      c.complex = 0;
      athin = isequalwithequalnans(a.inf,a.sup);
      bthin = isequalwithequalnans(b.inf,b.sup);
      if athin & bthin                    % R * R  with directed rounding
        setround(-1)
        c.inf = a.inf*b.inf;
        setround(1)
        c.sup = a.inf*b.inf;
      elseif athin
        setround(1)                       % R * IR
        bmid = b.inf + 0.5*(b.sup-b.inf);
        brad = bmid - b.inf;
        crad = abs(a.inf) * brad;
        c.inf = 0;                        % preserve order of definition of intval c
        c.sup = a.inf * bmid + crad;
        setround(-1)
        c.inf = a.inf * bmid - crad;
      elseif bthin
        setround(1)                       % IR * R
        amid = a.inf + 0.5*(a.sup-a.inf);
        arad = amid - a.inf;
        crad = arad * abs(b.inf);
        c.inf = 0;                        % preserve order of definition of intval c
        c.sup = amid * b.inf + crad;
        setround(-1)
        c.inf = amid * b.inf - crad;
      else                                  % IR * IR
        if INTLAB_INTVAL_IVMULT | (n==1)    % interval mtimes interval by outer product
          % index sets for first factor (also necessary for both factors containing proper zero intervals)
          index_ainf_neg = (a.inf<0);       % entries containing negative numbers
          index_asup_pos = (a.sup>0);       % entries containing positive numbers
          index_a_0 = (index_ainf_neg & index_asup_pos);   % proper zero intervals
          % index sets for second factor (also necessary for both factors containing proper zero intervals)
          index_binf_neg = (b.inf<0);       % entries containing negative numbers
          index_bsup_pos = (b.sup>0);       % entries containing positive numbers
          index_b_0 = (index_binf_neg & index_bsup_pos);   % proper zero intervals
          if ~any(any(index_a_0))           % first factor does not contain proper zero intervals
            % additional index sets for second factor
            index_binf_pos = (b.inf>0);     % positive entries 
            index_bsup_neg = (b.sup<0);     % negative entries
            if any(index_ainf_neg(:))       % there are non-positive intervals
              if issparse(a.inf)
                dd = a.sup;
                dd(a.inf<0) = -1;
                dd(a.sup>0) = 1;
                index = ( dd<0 );
                ainf = sparse([],[],[],m,n);
                ainf(index) = a.inf(index);
              else
                ainf = a.inf;
                ainf(index_asup_pos) = 0;   % only non-positive entries
              end
              if issparse(a.sup)
                [ia,ja,sa] = find(a.sup);
                index = ( sa<=0 );
                asup = sparse(ia(index),ja(index),sa(index),m,n);
              else
                asup = a.sup;
                asup(index_asup_pos) = 0;   % only non-positive entries
              end
              setround(-1)                  % lower bound, factor is b.sup
              if any(index_bsup_neg(:))
                if issparse(b.sup)
                  [ia,ja,sa] = find(b.sup);
                  index = ( sa<=0 );
                  bsup = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  bsup = b.sup;
                  bsup(index_bsup_pos) = 0;
                end
                c.inf = asup*bsup;
              else
                if issparse(asup)
                  c.inf = sparse([],[],[],m,n2);
                else
                  c.inf = zeros(m,n2);
                end
              end
              if any(index_bsup_pos(:))
                if issparse(b.sup)
                  [ia,ja,sa] = find(b.sup);
                  index = ( sa>=0 );
                  bsup = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  bsup = b.sup;
                  bsup(index_bsup_neg) = 0;
                end
                c.inf = c.inf + ainf*bsup;
              end
              setround(1)                 % upper bound, factor is b.inf
              if any(index_binf_neg(:))
                if issparse(b.inf)
                  [ia,ja,sa] = find(b.inf);
                  index = ( sa<=0 );
                  binf = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  binf = b.inf;
                  binf(index_binf_pos) = 0;
                end
                c.sup = ainf*binf;
              else
                if issparse(ainf)
                  c.sup = sparse([],[],[],m,n2);
                else
                  c.sup = zeros(m,n2);
                end
              end
              if any(index_binf_pos(:))
                if issparse(b.inf)
                  [ia,ja,sa] = find(b.inf);
                  index = ( sa>=0 );
                  binf = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  binf = b.inf;
                  binf(index_binf_neg) = 0;
                end
                c.sup = c.sup + asup*binf;
              end
            else
              if issparse(a.inf)
                c.inf = sparse([],[],[],m,n2);
              else                
                c.inf = zeros(m,n2);
              end
              c.sup = c.inf;
            end
            if any(index_asup_pos(:))     % there are non-positive intervals
              if issparse(a.inf)
                [ia,ja,sa] = find(a.inf);
                index = ( sa>=0 );
                ainf = sparse(ia(index),ja(index),sa(index),m,n);
              else
                ainf = a.inf;
                ainf(index_ainf_neg) = 0;   % only non-negative entries
              end
              if issparse(a.sup)
                dd = a.inf;
                dd(a.sup>0) = 1;
                dd(a.inf<0) = -1;
                index = ( dd>0 );
                asup = sparse([],[],[],m,n);
                asup(index) = a.sup(index);
              else
                asup = a.sup;
                asup(index_ainf_neg) = 0;   % only non-negative entries
              end
              setround(-1)                % lower bound
              if any(index_binf_neg(:))
                if issparse(b.inf)
                  [ia,ja,sa] = find(b.inf);
                  index = ( sa<=0 );
                  binf = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  binf = b.inf;
                  binf(index_binf_pos) = 0;
                end
                c.inf = c.inf + asup*binf;
              end
              if any(index_binf_pos(:))
                if issparse(b.inf)
                  [ia,ja,sa] = find(b.inf);
                  index = ( sa>=0 );
                  binf = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  binf = b.inf;
                  binf(index_binf_neg) = 0;
                end
                c.inf = c.inf + ainf*binf;
              end
              setround(1)                 % upper bound
              if any(index_bsup_neg(:))
                if issparse(b.sup)
                  [ia,ja,sa] = find(b.sup);
                  index = ( sa<=0 );
                  bsup = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  bsup = b.sup;
                  bsup(index_bsup_pos) = 0;
                end
                c.sup = c.sup + ainf*bsup;
              end
              if any(index_bsup_pos(:))
                if issparse(b.sup)
                  [ia,ja,sa] = find(b.sup);
                  index = ( sa>=0 );
                  bsup = sparse(ia(index),ja(index),sa(index),n,n2);
                else
                  bsup = b.sup;
                  bsup(index_bsup_neg) = 0;
                end
                c.sup = c.sup + asup*bsup;
              end
            end
            if ~issparse(c)               % take care of NaN
              index = isnan(a.inf) | isnan(a.sup);
              if any(index(:))
                c.inf(any(index,2),:) = NaN;
                c.sup(any(index,2),:) = NaN;
              end
              index = isnan(b.inf) | isnan(b.sup);
              if any(index(:))
                c.inf(:,any(index,1)) = NaN;
                c.sup(:,any(index,1)) = NaN;
              end
            end
          elseif ~any(any(index_b_0))     % second factor does not contain proper zero intervals
            % zero_a=1, i.e. a contains proper zero intervals
            % additional index sets for second factor
            index_ainf_pos = (a.inf>0);       % positive entries 
            index_asup_neg = (a.sup<0);       % negative entries
            if any(index_binf_neg(:))     % there are non-positive intervals
              if issparse(b.inf)
                dd = b.sup;
                dd(b.inf<0) = -1;
                dd(b.sup>0) = 1;
                index = ( dd<0 );
                binf = sparse([],[],[],n,n2);
                binf(index) = b.inf(index);
              else
                binf = b.inf;
                binf(index_bsup_pos) = 0;   % only non-positive entries
              end
              if issparse(b.sup)
                [ia,ja,sa] = find(b.sup);
                index = ( sa<=0 );
                bsup = sparse(ia(index),ja(index),sa(index),n,n2);
              else
                bsup = b.sup;
                bsup(index_bsup_pos) = 0;   % only non-positive entries
              end
              setround(-1)                % lower bound, factor is a.sup
              if any(index_asup_neg(:))
                if issparse(a.sup)
                  [ia,ja,sa] = find(a.sup);
                  index = ( sa<=0 );
                  asup = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  asup = a.sup;
                  asup(index_asup_pos) = 0;
                end
                c.inf = asup*bsup;
              else
                if issparse(bsup)
                  c.inf = sparse([],[],[],m,n2);
                else
                  c.inf = zeros(m,n2);
                end
              end
              if any(index_asup_pos(:))
                if issparse(a.sup)
                  [ia,ja,sa] = find(a.sup);
                  index = ( sa>=0 );
                  asup = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  asup = a.sup;
                  asup(index_asup_neg) = 0;
                end
                c.inf = c.inf + asup*binf;
              end
              setround(1)                 % upper bound, factor is a.inf
              if any(index_ainf_neg(:))
                if issparse(a.inf)
                  [ia,ja,sa] = find(a.inf);
                  index = ( sa<=0 );
                  ainf = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  ainf = a.inf;
                  ainf(index_ainf_pos) = 0;
                end
                c.sup = ainf*binf;
              else
                if issparse(binf)
                  c.sup = sparse([],[],[],m,n2);
                else
                  c.sup = zeros(m,n2);
                end
              end
              if any(index_ainf_pos(:))
                if issparse(a.inf)
                  [ia,ja,sa] = find(a.inf);
                  index = ( sa>=0 );
                  ainf = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  ainf = a.inf;
                  ainf(index_ainf_neg) = 0;
                end
                c.sup = c.sup + ainf*bsup;
              end
            else
              if issparse(a.inf)
                c.inf = sparse([],[],[],m,n2);
              else                
                c.inf = zeros(m,n2);
              end
              c.sup = c.inf;
            end
            if any(index_bsup_pos(:))     % there are non-positive intervals
              if issparse(b.inf)
                [ia,ja,sa] = find(b.inf);
                index = ( sa>=0 );
                binf = sparse(ia(index),ja(index),sa(index),n,n2);
              else
                binf = b.inf;
                binf(index_binf_neg) = 0;   % only non-negative entries
              end
              if issparse(b.sup)
                dd = b.inf;
                dd(b.sup>0) = 1;
                dd(b.inf<0) = -1;
                index = ( dd>0 );
                bsup = sparse([],[],[],n,n2);
                bsup(index) = b.sup(index);
              else
                bsup = b.sup;
                bsup(index_binf_neg) = 0;   % only non-negative entries
              end
              setround(-1)                % lower bound, factor is a.inf
              if any(index_ainf_neg(:))
                if issparse(a.inf)
                  [ia,ja,sa] = find(a.inf);
                  index = ( sa<=0 );
                  ainf = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  ainf = a.inf;
                  ainf(index_ainf_pos) = 0;
                end
                c.inf = c.inf + ainf*bsup;
              end
              if any(index_ainf_pos(:))
                if issparse(a.inf)
                  [ia,ja,sa] = find(a.inf);
                  index = ( sa>=0 );
                  ainf = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  ainf = a.inf;
                  ainf(index_ainf_neg) = 0;
                end
                c.inf = c.inf + ainf*binf;
              end
              setround(1)                 % upper bound, factor is a.sup
              if any(index_asup_neg(:))
                if issparse(a.sup)
                  [ia,ja,sa] = find(a.sup);
                  index = ( sa<=0 );
                  asup = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  asup = a.sup;
                  asup(index_asup_pos) = 0;
                end
                c.sup = c.sup + asup*binf;
              end
              if any(index_asup_pos(:))
                if issparse(a.sup)
                  [ia,ja,sa] = find(a.sup);
                  index = ( sa>=0 );
                  asup = sparse(ia(index),ja(index),sa(index),m,n);
                else
                  asup = a.sup;
                  asup(index_asup_neg) = 0;
                end
                c.sup = c.sup + asup*bsup;
              end
            end
            if ~issparse(c)               % take care of NaN
              index = isnan(a.inf) | isnan(a.sup);
              if any(index(:))
                c.inf(any(index,2),:) = NaN;
                c.sup(any(index,2),:) = NaN;
              end
              index = isnan(b.inf) | isnan(b.sup);
              if any(index(:))
                c.inf(:,any(index,1)) = NaN;
                c.sup(:,any(index,1)) = NaN;
              end
            end
          else                            % both factors contain proper zero intervals
            % index sets for first factor already computed
            % index_ainf_neg = ( a.inf < 0);  
            % index_asup_pos = ( a.sup > 0);
            % index_a_0 = (index_ainf_neg & index_asup_pos);
            % submatrices of first factor
            if issparse(a.inf)
              % avoid matrix(index)=0, very slow for sparse matrices with not so few elements
              % ainf_neg = a.inf; ainf_neg(index_asup_pos) = 0;         
              dd = a.sup;
              dd(a.inf<0) = -1;
              dd(a.sup>0) = 1;
              index = ( dd<0 );
              ainf_neg = sparse([],[],[],m,n);
              ainf_neg(index) = a.inf(index);
              % ainf_pos = a.inf; ainf_pos(index_ainf_neg) = 0;
              [ia,ja,sa] = find(a.inf);
              index = ( sa>=0 );
              ainf_pos = sparse(ia(index),ja(index),sa(index),m,n);
              % asup_neg = a.sup; asup_neg(index_asup_pos) = 0;
              [ia,ja,sa] = find(a.sup);
              index = ( sa<=0 );
              asup_neg = sparse(ia(index),ja(index),sa(index),m,n);
              % asup_pos = a.sup; asup_pos(index_ainf_neg) = 0;              
              dd = a.inf;
              dd(a.sup>0) = 1;
              dd(a.inf<0) = -1;
              index = ( dd>0 );
              asup_pos = sparse([],[],[],m,n);
              asup_pos(index) = a.sup(index);
              % prepare for ainf_0 and asup_0         
              ainf_0 = sparse([],[],[],m,n,0);
            else
              ainf_neg = a.inf; ainf_neg(index_asup_pos) = 0;         
              ainf_pos = a.inf; ainf_pos(index_ainf_neg) = 0;
              asup_neg = a.sup; asup_neg(index_asup_pos) = 0;
              asup_pos = a.sup; asup_pos(index_ainf_neg) = 0;
              % prepare for ainf_0 and asup_0   
              ainf_0 = zeros(m,n);
            end
            asup_0 = ainf_0;
            ainf_0(index_a_0) = a.inf(index_a_0); 
            asup_0(index_a_0) = a.sup(index_a_0);
            % index sets for second factor already computed
            % index_binf_neg = ( b.inf < 0); 
            % index_bsup_pos = ( b.sup > 0);
            % index_b_0 = (index_binf_neg & index_bsup_pos);
            % submatrices of second factor  
            if issparse(b.inf)
              % avoid matrix(index)=0, very slow for sparse matrices with not so few elements
              % binf_neg = b.inf; binf_neg(index_bsup_pos) = 0;
              dd = b.sup;
              dd(b.inf<0) = -1;
              dd(b.sup>0) = 1;
              index = ( dd<0 );
              binf_neg = sparse([],[],[],n,n2);
              binf_neg(index) = b.inf(index);
              % binf_pos = b.inf; binf_pos(index_binf_neg) = 0;
              [ia,ja,sa] = find(b.inf);
              index = ( sa>=0 );
              binf_pos = sparse(ia(index),ja(index),sa(index),n,n2);
              % bsup_neg = b.sup; bsup_neg(index_bsup_pos) = 0;
              [ia,ja,sa] = find(b.sup);
              index = ( sa<=0 );
              bsup_neg = sparse(ia(index),ja(index),sa(index),n,n2);
              % bsup_pos = b.sup; bsup_pos(index_binf_neg) = 0;
              dd = b.inf;
              dd(b.sup>0) = 1;
              dd(b.inf<0) = -1;
              index = ( dd>0 );
              bsup_pos = sparse([],[],[],n,n2);
              bsup_pos(index) = b.sup(index);
              % prepare for binf_0 and bsup_0         
              binf_0 = sparse([],[],[],n,n2,0);
            else
              binf_neg = b.inf; binf_neg(index_bsup_pos) = 0;
              binf_pos = b.inf; binf_pos(index_binf_neg) = 0;
              bsup_neg = b.sup; bsup_neg(index_bsup_pos) = 0;
              bsup_pos = b.sup; bsup_pos(index_binf_neg) = 0;
              % prepare for binf_0 and bsup_0   
              binf_0 = zeros(n,n2);
            end
            bsup_0 = binf_0;
            binf_0(index_b_0) = b.inf(index_b_0); 
            bsup_0(index_b_0) = b.sup(index_b_0);
            % check existence of both proper zero intervals
            index0 = find( any(index_a_0,1) & any(index_b_0',1) );
            % lower bound            
            setround(-1)
            % everything except proper zero intervals
            c.inf = asup_neg*bsup_neg + (ainf_neg+ainf_0)*bsup_pos + ainf_neg*bsup_0 + ...
              (asup_0+asup_pos)*binf_neg + asup_pos*binf_0 + ainf_pos*binf_pos;
            % lower bound of product of proper zero intervals
            if issparse(c.inf)              % sparse product
              delta = sparse(0);            % sparse initialization
              len0 = length(index0);
              % faster summation by unrolled loop adding sparse matrices with fewer elements
              N = 6;                        % terms in unrolled loop
              for kk=1:N:(len0-N+1)
                k = index0(kk);
                k1 = index0(kk+1);
                k2 = index0(kk+2);
                k3 = index0(kk+3);
                k4 = index0(kk+4);
                k5 = index0(kk+5);
                delta = delta + ...
                  ( min( ainf_0(:,k)*bsup_0(k,:) , asup_0(:,k)*binf_0(k,:) ) + ...
                  min( ainf_0(:,k1)*bsup_0(k1,:) , asup_0(:,k1)*binf_0(k1,:) ) + ...
                  min( ainf_0(:,k2)*bsup_0(k2,:) , asup_0(:,k2)*binf_0(k2,:) ) + ...
                  min( ainf_0(:,k3)*bsup_0(k3,:) , asup_0(:,k3)*binf_0(k3,:) ) + ...
                  min( ainf_0(:,k4)*bsup_0(k4,:) , asup_0(:,k4)*binf_0(k4,:) ) + ...
                  min( ainf_0(:,k5)*bsup_0(k5,:) , asup_0(:,k5)*binf_0(k5,:) ) );
              end
              % take care of possibly remaining terms
              for kk=(N*floor(len0/N)+1):len0
                k = index0(kk);
                delta = delta + min( ainf_0(:,k)*bsup_0(k,:) , asup_0(:,k)*binf_0(k,:) );
              end
              c.inf = c.inf + delta;
            else                            % full product
              for k=index0
                c.inf = c.inf + min( ainf_0(:,k)*bsup_0(k,:) , asup_0(:,k)*binf_0(k,:) );
              end
            end
            setround(1)
            % everything except proper zero intervals
            c.sup = asup_neg*binf_pos + (ainf_neg+ainf_0)*binf_neg + ainf_neg*binf_0 + ...
              (asup_0+asup_pos)*bsup_pos + asup_pos*bsup_0 + ainf_pos*bsup_neg;
            % upper bound of product of proper zero intervals
            if issparse(c.sup)              % sparse product
              delta = sparse(0);            % sparse initialization
              % faster summation by unrolled loop adding sparse matrices with fewer elements
              for kk=1:N:(len0-N+1)
                k = index0(kk);
                k1 = index0(kk+1);
                k2 = index0(kk+2);
                k3 = index0(kk+3);
                k4 = index0(kk+4);
                k5 = index0(kk+5);
                delta = delta + ...
                  ( max( ainf_0(:,k)*binf_0(k,:) , asup_0(:,k)*bsup_0(k,:) ) + ...
                  max( ainf_0(:,k1)*binf_0(k1,:) , asup_0(:,k1)*bsup_0(k1,:) ) + ...
                  max( ainf_0(:,k2)*binf_0(k2,:) , asup_0(:,k2)*bsup_0(k2,:) ) + ...
                  max( ainf_0(:,k3)*binf_0(k3,:) , asup_0(:,k3)*bsup_0(k3,:) ) + ...
                  max( ainf_0(:,k4)*binf_0(k4,:) , asup_0(:,k4)*bsup_0(k4,:) ) + ...
                  max( ainf_0(:,k5)*binf_0(k5,:) , asup_0(:,k5)*bsup_0(k5,:) ) );
              end
              % take care of possibly remaining terms
              for kk=(N*floor(len0/N)+1):len0
                k = index0(kk);
                delta = delta + max( ainf_0(:,k)*binf_0(k,:) , asup_0(:,k)*bsup_0(k,:) );
              end
              c.sup = c.sup + delta;
            else                            % full product
              for k=index0
                c.sup = c.sup + max( ainf_0(:,k)*binf_0(k,:) , asup_0(:,k)*bsup_0(k,:) );
              end
            end
          end
        else
          setround(1)                       % fast interval mtimes interval thru
          amid = a.inf + 0.5*(a.sup-a.inf);   % mid/rad arithmetic
          arad = amid - a.inf;
          bmid = b.inf + 0.5*(b.sup-b.inf);
          brad = bmid - b.inf;
          crad = arad * ( abs(bmid) + brad ) + abs(amid) * brad ;
          c.inf = 0;                    % preserve order of definition of intval c
          c.sup = amid*bmid + crad;
          setround(-1)
          c.inf = amid*bmid - crad;
        end
      end
      c.mid = [];
      c.rad = [];
    end
  end
  
  INTLAB_SPARSE_BUG = getappdata(0,'INTLAB_SPARSE_BUG');
  if INTLAB_SPARSE_BUG
    % take care of Matlab sparse NaN bug: sparse factor a
    if issparse(a)                            % simplified: factor inf always produces NaN
      index = any(isnan(b),1) | any(isinf(b),1);
      if any(any(index))
        if c.complex
          c.mid(:,index) = NaN;
          c.rad(:,index) = NaN;
        else
          c.inf(:,index) = NaN;
          c.sup(:,index) = NaN;
        end
      end
    end
    % take care of Matlab sparse NaN bug: sparse factor b
    if issparse(b)                            % simplified: factor inf always produces NaN
      index = any(isnan(a),2) | any(isinf(a),2);
      if any(any(index))
        if c.complex
          c.mid(index,:) = NaN;
          c.rad(index,:) = NaN;
        else
          c.inf(index,:) = NaN;
          c.sup(index,:) = NaN;
        end
      end
    end
  end
  
  % non-interval output for performance improvement for hessians
  if nargin==2
    c = class(c,'intval');
  end
  
  setround(rndold)
