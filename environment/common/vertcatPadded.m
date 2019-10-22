function [v] = vertcatPadded(v,u)

[v] = concatinate_padding(v,u);

end

function [v] = concatinate_padding(v,nl)

if iscell(nl)
    blk = {[]};
end

if ischar(nl)
    blk = ' ';
end

if size(v,2) > size(nl,2)
    % Get the horizontal buffer difference for vertical concatination
    bufferString = repmat(blk,[size(nl,1) abs(size(v,2) - size(nl,2))]);   % Get the difference in horezontal width
    % Expand the new text to match paragraph
    nl = [nl,bufferString];
else
    % Get the horizontal buffer difference for vertical concatination
    bufferString = repmat(blk,[size(v,1) abs(size(v,2) - size(nl,2))]);   % Get the difference in horezontal width
    % Expand the paragraph to match new text
    v = [v,bufferString];
end
% Concatinate the text to the end of the paragraph
v = [v;nl];
end
