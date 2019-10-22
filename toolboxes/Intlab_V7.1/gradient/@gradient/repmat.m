function a = repmat(a,m,n)
%REPMAT       Implements  repmat(a)  for gradients
%
%Functionality as in Matlab.
%

% written  05/03/09     S.M. Rump 
%

%Code adapted from Matlab R2008a
%

if nargin < 2
    error('repmat requires at least 2 inputs')
end

if nargin == 2
    if prod(size(m))==1
        siz = [m m];
    else
        siz = m;
    end
else
    siz = [m n];
end

if prod(size(a))==1
    %VVVV  a = a(ones(siz));
    s.type = '()'; s.subs = {ones(siz)}; a = subsref(a,s);
    %AAAA  Matlab bug fix
elseif ndims(a) == 2 & length(siz) == 2
    [m,n] = size(a);
    
    if (m == 1 & siz(2) == 1)
        %VVVV  a = a(ones(siz(1), 1), :);
        s.type = '()'; s.subs = {ones(siz(1), 1), ':'}; a = subsref(a,s);
        %AAAA  Matlab bug fix
    elseif (n == 1 & siz(1) == 1)
        %VVVV  a = a(:, ones(siz(2), 1));
        s.type = '()'; s.subs = {':', ones(siz(2), 1)}; a = subsref(a,s);
        %AAAA  Matlab bug fix
    else
        mind = (1:m)';
        nind = (1:n)';
        mind = mind(:,ones(1,siz(1)));
        nind = nind(:,ones(1,siz(2)));
        %VVVV  a = a(mind,nind);
        s.type = '()'; s.subs = {mind,nind}; a = subsref(a,s);
        %AAAA  Matlab bug fix
    end
else
    asiz = size(a);
    asiz = [asiz ones(1,length(siz)-length(asiz))];
    siz = [siz ones(1,length(asiz)-length(siz))];
    for i=length(asiz):-1:1
        ind = (1:asiz(i))';
        subs{i} = ind(:,ones(1,siz(i)));
    end
    %VVVV  a = a(subs{:});
    s.type = '()'; s.subs = {subs{:}}; a = subsref(a,s);
    %AAAA  Matlab bug fix
end
