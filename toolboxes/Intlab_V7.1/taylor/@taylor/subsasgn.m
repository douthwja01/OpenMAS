function r = subsasgn(r,s,b)
%SUBSASGN     Implements subscripted assignments for Taylor
%
%  example  r(2,:) = b
%

% written  05/21/09     S.M. Rump
% modified 02/28/10     S.M. Rump  multiple index
% modified 08/26/12     S.M. Rump  global variables removed
% modified 09/26/12     S.M. Rump  index handling (thanks to Matthew Weinstein)
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  K1 = getappdata(0,'INTLAB_TAYLOR_ORDER') + 1;

  if length(s)>1
    error('multiple indexing for Taylor assignment not allowed')
  end

  if strcmp(s.type,'()')     % assignment r(i) = b

    rEmpty = isempty(r);    % assignment r(i) = b for empty r
    if isempty(b)
      % does not work in Matlab 5.3 for sparse r
      index = reshape(1:prod(r.size),r.size);
      index(s.subs{:}) = [];
      r.size = size(index);
      r.t = r.t(:,index);
      if rndold
        setround(rndold)
      end
      return
    end

    if ~isa(b,'taylor')
      b = taylor(b);
    end

    if ~rEmpty
      resultIsintval = isa(r.t,'intval');
      if ~resultIsintval & isa(b.t,'intval')
        r = intval(r);
        resultIsintval = 1;
      end
      sizeincreased = 0;
      rsizeold = r.size;
      if length(s.subs)==1               % single index
        if ~isequal(s.subs{1},':')       % not call r(:)=...
          if r.size(1)==1                % row vector
            M = s.subs{1};
            sizeincreased = ( M>r.size(2) );
            if sizeincreased
              r.size(2) = M;
            end
          else
            sizeincreased = ( s.subs{1} > prod(r.size) );
            if sizeincreased
              error('In an assignment A(I) = B a Taylor matrix cannot be resized')
            end
          end
          if sizeincreased
            srx = r.size;
            if length(srx)==2
              if all(prod(srx)~=srx)
                error('matrix cannot be resized by assignment a(I) = b')
              end
            else
              error('attempt to grow size of array along ambiguous dimension')
            end
          end
        end
      else                            % multiple index
        for i=1:length(s.subs)
          if ~isequal(s.subs{i},':')
            M = max(s.subs{i});
            if ( M > r.size(i) )
              sizeincreased = 1;
              r.size(i) = M;
            end
          end
        end
      end
      if sizeincreased                % size increased, adapt .t
        rt = r.t;
        value = ones(rsizeold);
        value( s.subs{:} ) = 0;
        r.t = zeros(K1,prod(size(value)));
        if resultIsintval
          r.t = intval(r.t);
        end
        r.t( : , value==1 ) = rt;
      end
    else                     % assignment r(i) = b for empty r
      resultIsintval = isa(b.t,'intval');
      for i=1:length(s.subs)
        r.size(i) = max(s.subs{i});
      end
      if length(r.size)==1
        r.size = [1 r.size];
      end
      N = prod(r.size);
      if N>10
        r.t = sparse([],[],[],K1,N);
      else
        r.t = zeros(K1,N);
      end
      index = reshape(1:prod(r.size),r.size);
      r.t(:,index(s.subs{:})) = 0;
      if resultIsintval
        r.t = intval(r.t);
      end
    end
    value = reshape(1:prod(r.size),r.size);
    index = value( s.subs{:} );
    if ( prod(b.size)==1 ) & ( length(index)~=1 )
      r.t(:,index) = repmat(b.t,1,length(index));
    else
      r.t(:,index) = b.t;
    end
    if rEmpty
      r = class(r,'taylor');
    end
  else
    error('invalid index reference for taylor')
  end
    
  if rndold
    setround(rndold)
  end
