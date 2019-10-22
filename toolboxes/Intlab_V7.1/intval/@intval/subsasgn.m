function c = subsasgn(c,s,b)
%SUBSASGN     Implements subscripted assignments for intervals
%
%   example   c(i,:) = b
%

% written  10/16/98     S.M. Rump
% modified 11/23/98     S.M. Rump  delete components by a(i,j) = []
% modified 07/01/99     S.M. Rump  multi-dimensional arrays, complex=0
% modified 09/02/00     S.M. Rump  real-complex assignment
% modified 09/28/01     S.M. Rump  empty left or right hand side
% modified 09/29/02     S.M. Rump  sparse parameters
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    remove check for 'double'
%                                    Matlab sparse bug
%                                    complex assignment a(n,n)=.. for huge n
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 10/13/12     S.M. Rump  INTLAB_INTVAL_STDFCTS
%

if length(s)>1
    error('multiple indexing for interval assignment not allowed')
end

if strcmp(s.type,'()')     % subarray assignment c(i,j) = b
    % on entry, l.h.s. c is either of type interval or, not defined.
    % The latter case is equivalent to isa(c,'intval')=0, because
    % always exist('c')=1.
    
    if iscell(s(1).subs) == 1 && iscell(s(1).subs{1}) == 1
    	s(1).subs = s(1).subs{1};
    end
    
    if ~isa(b,'intval')
        b = intval(b);
    end
    if ~isa(c,'intval')
        c.complex = b.complex;
        c.inf = [];
        c.sup = [];
        c.mid = [];
        c.rad = [];
        c = class(c,'intval');
    end
    INTLAB_STDFCTS_RCASSIGN = getappdata(0,'INTLAB_STDFCTS_RCASSIGN');
    if INTLAB_STDFCTS_RCASSIGN
        if c.complex & ~b.complex
            b = cintval(b);
        end
        if ~c.complex & b.complex
            if INTLAB_STDFCTS_RCASSIGN==1
                warning('**** Subscripted assignment  real(...) = complex')
            else
                error('**** Subscripted assignment  real(...) = complex')
            end
            c = cintval(c);
        end
    else
        if b.complex~=c.complex
            if b.complex
                c = cintval(c);
            else
                b = cintval(b);
            end
        end
    end
    if c.complex
        if isequal(c.rad,0)
            if issparse(c.mid)
                c.rad = sparse([],[],[],size(c.mid,1),size(c.mid,2));
            else
                c.rad = zeros(size(c.mid));
            end
        end
        if isempty(b.mid)
            c.mid(s.subs{:}) = [];
            c.rad(s.subs{:}) = [];
        else
            c.mid(s.subs{:}) = b.mid;
            c.rad(s.subs{:}) = b.rad;
        end
        if ~any(c.rad)
            c.rad = 0;
        end
    else
        if isempty(b.inf)
            c.inf(s.subs{:}) = [];
            c.sup(s.subs{:}) = [];
        else
            c.inf(s.subs{:}) = b.inf;
            c.sup(s.subs{:}) = b.sup;
        end
    end
elseif strcmp(s.type,'.')      % subfield access
    if strcmp(s.subs,'inf')
        error('for safety, explicit assignment of x.inf=y not allowed; use x=infsup(y,sup(x))')
    elseif strcmp(s.subs,'sup')
        error('for safety, explicit assignment of x.sup=y not allowed; use x=infsup(inf(x),y)')
    elseif strcmp(s.subs,'mid')
        error('for safety, explicit assignment of x.mid=y not allowed; use x=midrad(y,rad(x))')
    elseif strcmp(s.subs,'rad')
        error('for safety, explicit assignment of x.rad=y not allowed; use x=midrad(mid(x),y)')
    else
        error('invalid call of subsasgn')
    end
else
    error('invalid call of subsasgn')
end

% avoid Matlab 6.5f bug:
% a = sparse([],[],[],1,1); a = reshape(a,1,1); abs(a)
% produces  9.6721e-317  or similar number in underflow range
if c.complex
    if prod(size(c.mid))==1
        c.mid = full(c.mid);
        c.rad = full(c.rad);
    end
else
    if prod(size(c.inf))==1
        c.inf = full(c.inf);
        c.sup = full(c.sup);
    end
end
