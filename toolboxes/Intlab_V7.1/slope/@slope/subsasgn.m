function u = subsasgn(u,s,b)
%SUBSASGN     Implements subscripted assignments for slopes
%
%  example  u(2,:) = b
%

% written  12/06/98     S.M. Rump
% modified 10/15/99     S.M. Rump  assignment []
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 09/29/02     S.M. Rump  fix due to different behaviour of logical of Matlab 6.5
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    sparse input
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 09/26/12     S.M. Rump  index handling (thanks to Matthew Weinstein)
%

  if length(s)>1
    error('multiple indexing for slope assignment not allowed')
  end

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  if strcmp(s.type,'()')     % assignment u(i) = b

    uEmpty = isempty(u);     % assignment u(i) = b for empty u
    if isempty(b)
      % does not work in Matlab 5.3 for sparse r
      index = reshape(1:prod(u.size),u.size);
      index(s.subs{:}) = [];
      u.size = size(index);
      u.r = u.r(index,:);
      u.s = u.s(index,:);
      return
    end

    if ~isa(b,'slope')
      b = slope(b);
    end

    if ~uEmpty
      
      sizeincreased = 0;
      if length(s.subs)==1               % single index
        if ~isequal(s.subs{1},':')       % not call r(:)=...
          if u.size(1)==1                % row vector
            sizeincreased = ( s.subs{1}>u.size(2) );
          else
            sizeincreased = ( s.subs{1} > prod(u.size) );
          end
          if sizeincreased
            sux = u.size;
            if length(sux)==2
              if all(prod(sux)~=sux)
                error('matrix cannot be resized by assignment a(I) = b')
              end
            else
              error('attempt to grow size of array along ambiguous dimension')
            end
          end
        end
      else                            % multiple index
        for i=1:length(s.subs)
          if ( ~isequal(s.subs{i},':') ) & ( s.subs{i}~=1 )
            if ( i>length(u.size) ) 
              sizeincreased = 1;
            else
              sizeincreased = sizeincreased | ( s.subs{i} > u.size(i) );
            end
          end
        end
      end
      if sizeincreased                % size increased, adapt .r and .s
        ur = u.r;
        us = u.s;
        value = ones(u.size);
        value( s.subs{:} ) = 0;
        u.size = size(value);
        if issparse(ur)
          u.r = intval(sparse([],[],[],prod(u.size),INTLAB_SLOPE.NUMVAR+1));
          u.s = intval(sparse([],[],[],prod(u.size),INTLAB_SLOPE.NUMVAR));
        else
          u.r = intval(zeros(prod(u.size),INTLAB_SLOPE.NUMVAR+1));
          u.s = intval(zeros(prod(u.size),INTLAB_SLOPE.NUMVAR));
        end
        index = ( value(:)==1 );
        u.r(index,:) = ur;
        u.s(index,:) = us;
      end
      
      value = reshape(1:prod(u.size),u.size);
      index = value( s.subs{:} );
      if ( prod(b.size)==1 ) & ( length(index)~=1 )
        u.r(index,:) = repmat(b.r,length(index),1);
        u.s(index,:) = repmat(b.s,length(index),1);
      else
        u.r(index,:) = b.r;
        u.s(index,:) = b.s;
      end
    else                    % assignment u(i) = b for empty u
      index(s.subs{:}) = 1;
      u.size = size(index);
      if issparse(b.s)
        u.r = intval(sparse([],[],[],prod(size(index)),INTLAB_SLOPE.NUMVAR+1));;
        u.s = intval(sparse([],[],[],prod(size(index)),INTLAB_SLOPE.NUMVAR));;
      else
        u.r = intval(zeros(prod(size(index)),INTLAB_SLOPE.NUMVAR+1));
        u.s = intval(zeros(prod(size(index)),INTLAB_SLOPE.NUMVAR));
      end
      if prod(b.size)==1
        u.r(index==1,:) = b.r(ones(sum(index(:)),1),:);
        u.s(index==1,:) = b.s(ones(sum(index(:)),1),:);
      else
        u.r(index==1,:) = b.r;
        u.s(index==1,:) = b.s;
      end
      u = class(u,'slope');
    end

  else
    error('invalid index reference for slope')
  end
