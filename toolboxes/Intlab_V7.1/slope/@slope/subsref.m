function u = subsref(a,s)
%SUBSREF      Implements subscripted references for slopes
%

% written  12/06/98     S.M. Rump
% modifed  09/14/00     S.M. Rump  a.s(i,j) fixed
% modified 09/28/01     S.M. Rump  matrices and multi-dimensional arrays
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    access of .s of sparse arrays, .mid allowed
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/05/12     S.M. Rump  internal use
%

  INTLAB_SLOPE = getappdata(0,'INTLAB_SLOPE');

  while 1
    if ~isa(a,'slope')
      u = subsref(a,s(1));
    elseif strcmp(s(1).type,'()')        % index reference a(i)
      index = reshape(1:prod(a.size),a.size);
      index = index(s(1).subs{:});
      u.size = size(index);
      u.r = a.r(index(:),:);
      u.s = a.s(index(:),:);
      u = class(u,'slope');
    elseif strcmp(s(1).type,'.')         % index reference a.c, a.r or a.s
      if strcmp(s(1).subs,'c')
        u = reshape(a.r(:,1),a.size);
      elseif strcmp(s(1).subs,'r')
        u = reshape(a.r(:,size(a.r,2)),a.size);
      elseif strcmp(s(1).subs,'s')
        if length(a.size)==2 & ( a.size(2)==1 )
          u = reshape(a.s,[a.size(1) size(a.s,2)]);
        else
          if issparse(a.s) & ( prod(a.size(2:end))>1 )
            error('access of .s of sparse arrays with more than one column, see slopeinit')
          end
          u = reshape(a.s,[a.size size(a.s,2)]);
        end
      elseif strcmp(s(1).subs,'ss')
        u = a.s;
      elseif strcmp(s(1).subs,'size')
        u = a.size;
      elseif strcmp(s(1).subs,'mid')
        u = mid(a);
      else
        error('invalid subscript reference for slope')
      end
    else
      error('invalid index reference for slope')
    end
    if length(s)==1
      return
    end
    s = s(2:end);
    a = u;
  end
