function c = times(a,b,dummy)
%TIMES        Implements  a .* b  for intervals
%

% written  10/16/98     S.M. Rump
% modified 09/02/00     S.M. Rump  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    remove check for 'double'
%                                    extra argument for non-interval output (internal use only)
%                                    take care of Matlab sparse Inf/NaN bug
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%                                    sparse radius for complex data
% modified 05/23/06     S.M. Rump  sparse Inf/NaN bug corrected in Version 7.2+
% modified 12/03/06     S.M. Rump  Sparse Bug global flag (thanks to Arnold)
% modified 12/05/06     S.M. Rump  0*infsup(0,inf) (thanks to Arnold)
% modified 01/19/07     S.M. Rump  zero radius for huge arrays
% modified 05/22/07     S.M. Rump  check for sparse zero factor
% modified 10/21/07     S.M. Rump  Matlab bug for sparse(0)*inf
% modified 10/18/08     S.M. Rump  abss replaced by mag, improved performance
% modified 08/26/12     S.M. Rump  global variables removed
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
  end
  huge = ( max(prod(size(a)),prod(size(b)))>1e9 );

  if ~isa(a,'intval')                     % a is double
    if ~isreal(a) | b.complex             % complex case
      if ~b.complex
        setround(1)
        b.mid = b.inf + 0.5*(b.sup-b.inf);
        b.rad = b.mid - b.inf;
      end
      c.complex = 1;
      c.inf = [];
      c.sup = [];
      if isreal(a) | ~b.complex           % R .* IC  or  C .* IC
        setround(-1)
        c1 = a .* b.mid;
        setround(1)
        c2 = a .* b.mid;
      else                                % C .* IC
        setround(-1)
        c1 = real(a) .* real(b.mid) + (-imag(a)) .* imag(b.mid) + ...
             ( real(a) .* imag(b.mid) + imag(a) .* real(b.mid) ) * j;
        setround(1)
        c2 = real(a) .* real(b.mid) + (-imag(a)) .* imag(b.mid) + ...
             ( real(a) .* imag(b.mid) + imag(a) .* real(b.mid) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % R .* IC  or  C .* IC
      c.rad = abs(c.mid-c1);
      if ~isequal(b.rad,0)                % make sure c.rad remains sparse
        c.rad = c.rad + abs(a) .* b.rad;
      end
      % too rare, improvement doesn't pay
%       if huge
%         [cradi,cradj] = find(c.rad);
%         if ~any(find(cradi))              % take care of huge arrays
%           c.rad = 0;
%         end
%       else
%         if ~any(find(c.rad))
%           c.rad = 0;
%         end
%       end
    else                                  % real case  R .* IR
      c.complex = 0;
      if ( prod(size(a))==1 ) & isequal(a(1,1),0)    % careful with 0*inf and sparse a
        if isinf(b.inf) | isinf(b.sup)
          c.inf = repmat(NaN,size(b.inf));
        else
          if issparse(b.inf)
            [mb,nb] = size(b.inf);
            c.inf = sparse([],[],[],mb,nb);
          else
            c.inf = zeros(size(b.inf));
          end
        end
        c.sup = c.inf;
      else
        setround(-1)
        c.inf = min( a .* b.inf , a .* b.sup );
        setround(1)
        c.sup = max( a .* b.inf , a .* b.sup );
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
      c.complex = 1;
      c.inf = [];
      c.sup = [];
      if ~a.complex | isreal(b)           % IC .* R  or  IC .* C
        setround(-1)
        c1 = a.mid .* b;
        setround(1)
        c2 = a.mid .* b;
      else                                % IC .* C
        setround(-1)
        c1 = real(a.mid) .* real(b) + (-imag(a.mid)) .* imag(b) + ...
             ( real(a.mid) .* imag(b) + imag(a.mid) .* real(b) ) * j;
        setround(1)
        c2 = real(a.mid) .* real(b) + (-imag(a.mid)) .* imag(b) + ...
             ( real(a.mid) .* imag(b) + imag(a.mid) .* real(b) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % R .* IC  or  C .* IC
      c.rad = abs(c.mid-c1);
      if ~isequal(a.rad,0)                % make sure c.rad remains sparse
        c.rad = c.rad + a.rad .* abs(b) ;
      end
      % too rare, improvement doesn't pay
%       if huge
%         [cradi,cradj] = find(c.rad);
%         if ~any(find(cradi))              % take care of huge arrays
%           c.rad = 0;
%         end
%       else
%         if ~any(find(c.rad))
%           c.rad = 0;
%         end
%       end
    else                                  % real case  IR .* R
      c.complex = 0;
      setround(-1)
      c.inf = min( a.inf .* b , a.sup .* b );
      setround(1)
      c.sup = max( a.inf .* b , a.sup .* b );
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
        c1 = a.mid .* b.mid;
        setround(1)
        c2 = a.mid .* b.mid;
      else
        setround(-1)
        c1 = real(a.mid) .* real(b.mid) + (-imag(a.mid)) .* imag(b.mid) + ...
             ( real(a.mid) .* imag(b.mid) + imag(a.mid) .* real(b.mid) ) * j;
        setround(1)
        c2 = real(a.mid) .* real(b.mid) + (-imag(a.mid)) .* imag(b.mid) + ...
             ( real(a.mid) .* imag(b.mid) + imag(a.mid) .* real(b.mid) ) * j;
      end
      c.mid = c1 + 0.5*(c2-c1);           % IR .* IC  or  IC .* IC
      c.rad = abs(c.mid-c1);
      if ~isequal(a.rad,0)                % make sure c.rad remains sparse
        if ~isequal(b.rad,0)
          c.rad = c.rad + a.rad .* ( abs(b.mid) + b.rad );
        else
          c.rad = c.rad + a.rad .* abs(b.mid);
        end
      end
      if ~isequal(b.rad,0)
        c.rad = c.rad + abs(a.mid) .* b.rad;
      end
      % too rare, improvement doesn't pay
%       if huge
%         [cradi,cradj] = find(c.rad);
%         if ~any(find(cradi))              % take care of huge arrays
%           c.rad = 0;
%         end
%       else
%         if ~any(find(c.rad))
%           c.rad = 0;
%         end
%       end
    else                                  % real case  IR .* IR
      c.complex = 0;
      %setround(-1)
      c.inf = min( a.inf .* b.inf , a.inf .* b.sup );
      c.inf = min( c.inf , a.sup .* b.inf );
      c.inf = min( c.inf , a.sup .* b.sup );
      %setround(1)
      c.sup = max( a.inf .* b.inf , a.inf .* b.sup );
      c.sup = max( c.sup , a.sup .* b.inf );
      c.sup = max( c.sup , a.sup .* b.sup );
      c.mid = [];
      c.rad = [];
    end
  end
  
  INTLAB_SPARSE_BUG = getappdata(0,'INTLAB_SPARSE_BUG');
  if INTLAB_SPARSE_BUG
    % take care of Matlab sparse NaN bug
    if huge                                       % huge linear indices, for simplicity inf*x=NaN
      if c.complex
        [m,n] = size(c.mid);
      else
        [m,n] = size(c.inf);
      end
      if issparse(a)
        index = any( isnan(b) | isinf(b) );
        if any(index(:))                          % keep linear index small
          [indexi,indexj] = find(isnan(b));
          cNaN = sparse(indexi,indexj,NaN,m,n);
          [ia,ja,aa] = find(mag(a));
          [ib,jb,bb] = find(mag(b));
          ab = sparse([ia;ib],[ja;jb],[complex(aa(:),0);complex(0,bb(:))],m,n);
          [indexi,indexj,ab] = find(ab);
          index = isinf(imag(ab));
          ab = real(ab);
          ab(index & (ab==0)) = NaN;
          ab(~isnan(ab)) = 0;
          cNaN = cNaN + sparse(indexi,indexj,ab,m,n);
          if c.complex
            c.mid = c.mid + cNaN;
            if isequal(c.rad,0)
              c.rad = cNaN;                       % scalar + sparse = full also for scalar=0
            else
              c.rad = c.rad + cNaN;
            end
          else
            c.inf = c.inf + cNaN;
            c.sup = c.sup + cNaN;
          end
        end
      end
      if issparse(b)
        index = any( isnan(a) | isinf(a) );
        if any(index(:))                          % keep linear index small
          [indexi,indexj] = find(isnan(a));
          cNaN = sparse(indexi,indexj,NaN,m,n);
          [ia,ja,aa] = find(mag(a));
          [ib,jb,bb] = find(mag(b));
          ab = sparse([ia;ib],[ja;jb],[complex(0,aa(:));complex(bb(:),0)],m,n);
          [indexi,indexj,ab] = find(ab);
          index = isinf(imag(ab));
          ab = real(ab);
          ab(index & (ab==0)) = NaN;
          ab(~isnan(ab)) = 0;
          cNaN = cNaN + sparse(indexi,indexj,ab,m,n);
          if c.complex
            c.mid = c.mid + cNaN;
            if isequal(c.rad,0)
              c.rad = cNaN;                       % scalar + sparse = full also for scalar=0
            else
              c.rad = c.rad + cNaN;
            end
          else
            c.inf = c.inf + cNaN;
            c.sup = c.sup + cNaN;
          end
        end
      end
    else
      Index = logical(0);
      if issparse(a) & ( prod(size(a))~=1 )     % sparse scalar is handled correctly by Matlab
        Index = isnan(b);
        if any(Index(:))
          if prod(size(b))==1
            Index = find( a==0 );					      % may be almost full
          end
        end
        index = find(isinf(b));
        if any(index)
          if prod(size(b))==1                   % a is scalar
            index = find( a==0 );               % may be almost full
          else                                  % a and b not scalar, same size
            if isa(a,'intval')
              if a.complex
                if isequal(a.rad,0)
                  index( a.mid(index)~=0 ) = [];
                else
                  index( ( a.mid(index)~=0 ) | ( a.rad(index)~=0 ) ) = [];
                end
              else
                index( ( a.inf(index)~=0 ) | ( a.sup(index)~=0 ) ) = [];
              end
            else
              index( a(index)~=0 ) = [];
            end
          end
          Index(index) = 1;
        end
      end
      if issparse(b) & ( prod(size(b))~=1 )     % sparse scalar is handled correctly by Matlab
        Index = Index | isnan(a);
        if any(Index(:))
          if prod(size(a))==1
            Index = find( b==0 );					      % may be almost full
          end
        end
        index = find(isinf(a));
        if any(index)
          if prod(size(a))==1                   % a is scalar
            index = find( b==0 );               % may be almost full
          else                                  % a and b not scalar, same size
            if isa(b,'intval')
              if b.complex
                if isequal(b.rad,0)
                  index( b.mid(index)==0 ) = [];
                else
                  index( ( b.mid(index)~=0 ) | ( b.rad(index)~=0 ) ) = [];
                end
              else
                index( ( b.inf(index)~=0 ) | ( b.sup(index)~=0 ) ) = [];
              end
            else
              index( b(index)~=0 ) = [];
            end
          end
          Index(index) = 1;
        end
      end
      if any(Index(:))
        if c.complex
          c.mid(Index) = NaN;
          if isequal(c.rad,0)
            [m,n] = size(c.mid);
            c.rad = sparse([],[],[],m,n);
          end
          c.rad(Index) = NaN;
        else
          c.inf(Index) = NaN;
          c.sup(Index) = NaN;
        end
      end
    end
  end
  
  % non-interval output for performance improvement for hessians
  if nargin==2
    c = class(c,'intval');
  end
  
  setround(rndold)
