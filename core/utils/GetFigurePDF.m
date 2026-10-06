%% Convert .fig to .pdf (GetFigurePDF.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isSuccessful,message] = GetFigurePDF(figureHandle,figurePath)

% Input sanity check
switch nargin
    case 0
        figureHandle = gcf;
        figurePath = [cd,'\',figureHandle.Name];
    case 1
        figurePath = [cd,'\',figureHandle.Name];
    otherwise
        % Do nothing
end
% Assertions
assert(ischar(figurePath),'Path must be a string');

% Attempt to push the figure to pdf
try
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(figureHandle,figurePath,'-dpdf','-r0','-bestfit');
    % Report
    message = ''; isSuccessful = true;
catch publishError
    % Report
    message = publishError.message;
    isSuccessful = false;
    % Rethrow for clarity
    warning(message);
end

end
