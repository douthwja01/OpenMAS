function c = subsref(a,s)
%SUBSREF      Implements subscripted references for intervals
%
%   example   c = a(:,3:5)
%

% written  10/16/98     S.M. Rump
% modified 11/23/98     S.M. Rump  improved speed
% modified 09/14/00     S.M. Rump  a.inf(i) fixed
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    Matlab sparse bug
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  while 1
    if ~isa(a,'intval')               % for a.inf(i) etc.
      c = subsref(a,s(1));
    elseif strcmp(s(1).type,'()')     % index reference a(i,j)
      c.complex = a.complex;
      if a.complex                    % a is complex
        c.inf = [];
        c.sup = [];
        c.mid = a.mid(s(1).subs{:});
        % careful: just 0 is not sparse and may cause tremendous memory need
        if isequal(a.rad,0)
          c.rad = a.rad;
        else
          c.rad = a.rad(s(1).subs{:});
        end
        % avoid Matlab 6.5f bug: 
        % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
        % produces  9.6721e-317  or similar number in underflow range
        if prod(size(c.mid))==1
          c.mid = full(c.mid);
        end
        if prod(size(c.rad))==1
          c.rad = full(c.rad);
        end
      else                          % a is real
        c.inf = a.inf(s(1).subs{:});
        c.sup = a.sup(s(1).subs{:});
        % avoid Matlab 6.5f bug: 
        % a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
        % produces  9.6721e-317  or similar number in underflow range
        if prod(size(c.inf))==1
          c.inf = full(c.inf);
        end
        if prod(size(c.sup))==1
          c.sup = full(c.sup);
        end
        c.mid = [];
        c.rad = [];
      end
      c = class(c,'intval');
    elseif strcmp(s(1).type,'.')     % subfield access
      if     strcmp(s(1).subs,'inf'), if a.complex, c=inf(a); else c=a.inf; end
      elseif strcmp(s(1).subs,'sup'), if a.complex, c=sup(a); else c=a.sup; end
      elseif strcmp(s(1).subs,'mid'), if a.complex, c=a.mid; else c=mid(a); end
      elseif strcmp(s(1).subs,'rad'), if a.complex, c=a.rad; else c=rad(a); end
      else
        error('invalid subscript reference for intval')
      end
    else
      error('invalid index reference for intval')
    end
    if length(s)==1
      return
    end
    s = s(2:end);
    a = c;
  end
