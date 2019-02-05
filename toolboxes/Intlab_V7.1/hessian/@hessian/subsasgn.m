function r = subsasgn(r,s,b)
%SUBSASGN     Implements subscripted assignments for hessians
%
%  example  r(2,:) = b
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 05/09/07     S.M. Rump  assignment r(:)=...
% modified 08/26/12     S.M. Rump  global variables removed
% modified 09/26/12     S.M. Rump  index handling (thanks to Matthew Weinstein)
%

  N = getappdata(0,'INTLAB_HESSIAN_NUMVAR');

  if length(s)>1
    error('multiple indexing for hessian assignment not allowed')
  end

  if strcmp(s.type,'()')     % assignment r(i) = b

    rEmpty = isempty(r);     % assignment r(i) = b for empty r
    if isempty(b)
      % does not work in Matlab 5.3 for sparse r
      r.x(s.subs{:}) = [];
      value = zeros(size(r.x));
      value(s.subs{:}) = 1;
      index = find(value);
      r.dx( : , index ) = [];
      r.hx( : , index ) = [];
      return
    end

    if ~isa(b,'hessian')
      b = hessian(b);
    end

    if ~rEmpty
      resultIsintval = isa(r.x,'intval');
      if ~resultIsintval & isa(b.x,'intval')
        r = intval(r);
        resultIsintval = 1;
      end
      sizeincreased = 0;
      if length(s.subs)==1              % single index
        if ~isequal(s.subs{1},':')      % not r(:)=...
          if size(r.x,1)==1             % row vector
            sizeincreased = ( s.subs{1}>size(r.x,2) );
          else
            sizeincreased = ( s.subs{1} > prod(size(r.x)) );
          end
          if sizeincreased
            srx = size(r.x);
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
            sizeincreased = sizeincreased | ( s.subs{i} > size(r.x,i) );
          end
        end
      end
      if sizeincreased                % size increased, adapt .x and .dx and .hx
        rx = r.x;
        rdx = r.dx;
        rhx = r.hx;
        value = ones(size(r.x));
        value( s.subs{:} ) = 0;
        r.x = zeros(size(value));
        INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
        if N<INTLAB_HESSIAN_SPARSE
          r.dx = zeros(N,prod(size(value)));
          r.hx = zeros(N*N,prod(size(value)));
        else
          r.dx = sparse([],[],[],N,prod(size(value)));
          r.hx = sparse([],[],[],N*N,prod(size(value)));
      end
        if resultIsintval
          r.x = intval(r.x);
          r.dx = intval(r.dx);
          r.hx = intval(r.hx);
        end
        index = find(value);
        r.x( index ) = rx;
        r.dx( : , index ) = rdx;
        r.hx( : , index ) = rhx;
      end
    else                     % assignment r(i) = b for empty r
      resultIsintval = isa(b.x,'intval');
      r.x(s.subs{:}) = 0;
      INTLAB_HESSIAN_SPARSE = getappdata(0,'INTLAB_HESSIAN_SPARSE');
      if N<INTLAB_HESSIAN_SPARSE
        r.dx = zeros(N,prod(size(r.x)));
        r.hx = zeros(N*N,prod(size(r.x)));
      else
        r.dx = sparse([],[],[],N,prod(size(r.x)));
        r.hx = sparse([],[],[],N*N,prod(size(r.x)));
      end
      if resultIsintval
        r.x = intval(r.x);
        r.dx = intval(r.dx);
        r.hx = intval(r.hx);
      end
    end
    r.x( s.subs{:} ) = b.x;
    value = reshape(1:prod(size(r.x)),size(r.x));
    index = value( s.subs{:} );
    if ( prod(size(b.x))==1 ) & ( length(index)~=1 )
      r.dx(:,index) = b.dx(:,ones(1,length(index)));
      r.hx(:,index) = b.hx(:,ones(1,length(index)));
    else
      r.dx(:,index) = b.dx;
      r.hx(:,index) = b.hx;
    end
    if rEmpty
      r = class(r,'hessian');
    end
  else
    error('invalid index reference for hessian')
  end

  % avoid Matlab 6.5f bug: 
  % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
  % produces  9.6721e-317  or similar number in underflow range
  if prod(size(r.x))==1
    r.x = full(r.x);
  end
  if prod(size(r.dx))==1
    r.dx = full(r.dx);
  end
  if prod(size(r.hx))==1
    r.hx = full(r.hx);
  end
  