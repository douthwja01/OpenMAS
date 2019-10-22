%% TABBED FIGURE GENERATOR (GetTabbedFigure.m) $%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = GetTabbedFigure(figHandles,groupName)

% Creates a tabbed window with a figure contained in each tab.
% 
% No tab workarounds necessary - MATLAB's built-in but undocumented tab
% features are used for smooth graphics handling and a clean interface.
%
% Input parameters:
%   figHandles: an array of handles to currently valid, visible figures.
%   These figures can be GUIDE based.
%
% Output parameters:
%   (optional): The handle to the tabbed figure is output if desired.
%
% Example: 
%     f1 = figure('Name','Sin Wave');
%     plot(sin(1:100));
%     f2 = figure('Name','Random Points');
%     plot(magic(5))
%     f3 = msgbox('A message box');
%     tabbedFig = figs2tabs([f1,f2,f3])
% 
% Warning:
%    This code heavily relies on undocumented Matlab functionality.
%    It has only been tested on Matlab 2012a. Use at your own risk.
%
% Known limitations: 
%   *uimenu's are not preserved
%   *figure-wide callback functions are not preserved, i.e. KeyPressFcn,
%   CloseRequestFcn.
%
% Known issues: 
%   *Not compatible with older versions of Matlab, though it could be with
%   some slight modifications to this function. Compatibility fixes to this
%   function are welcome from those with older versions.
%
% See also uitab, uitabgroup, setappdata, getappdata, guidata

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Chad Smith 

warning('off','MATLAB:uitabgroup:OldVersion')
if nargin < 1
    error('An array of figure handles must be input')
end
if nargin < 2 || ~exist('groupName','var')
    groupName = 'MATLAB Tabbed GUI';            % Apply default name
end
% create figure and tab group
[tabbedFig,tabGroupH] = init_figure(groupName);

% create tab for each figure
for tabNum=1:length(figHandles)
   figHandle = figHandles(tabNum);
   add_tab(tabbedFig, tabGroupH, figHandle, tabNum);
end

%set the first tab as active and update to ensure the right guidata, windowsize, etc.
initialize_tab_group(tabbedFig, tabGroupH);

%set output
if nargout == 1
    varargout{1} = tabbedFig;
end
end

function [tabbedFig,tabGroupH] = init_figure(figureName)
%create tabbed figure
tabbedFig = figure('Name',figureName, ... %title bar text
    'Tag','tabbedWindow',... 
    'NumberTitle','off', ... %hide figure number in title
    'IntegerHandle','off',... %use number like 360.0027 instead of 1,2, or 3
    'Resize','on'); %allow user to resize, TODO make contents normalized to allow for proportional resizing
    %'Menubar','none',... %dont have file menu, etc.

%create a tab group
tabGroupH = uitabgroup;
% set(tabGroupH,'SelectionChangeCallback',@update_guidata_and_resize) 
set(tabGroupH,'SelectionChangedFcn',@update_guidata_and_resize) 
drawnow
end

function initialize_tab_group(tabbedFig, tabGroupH)
%set the first tab as active and update position
curTabNum = 1;
% get(tabGroupH)
% pause;
% set(tabGroupH,'SelectedIndex',curTabNum)
%set(tabGroupH,'SelectedTab',curTabNum)

guiAndTabInfo =  getappdata(tabbedFig,'guiAndTabInfo');
newPos = calcNewPos(curTabNum, guiAndTabInfo);
set(tabbedFig,'Position',newPos)
CenterWindow(tabbedFig)
end

function add_tab(tabbedFig, tabGroupH, figHandle, tabNum)
%get all children a standalone figure
allChildren = get(figHandle,'Children');

%isolate type "uimenu"
%determine types of children

try
    types = get(allChildren,'Type');
    types = confirm_cell(types);
    uiMenuIndxsBool = cregexp(types,'uimenu');
    
    %add all children except those of type "uimenu"
    validChildren = allChildren(~uiMenuIndxsBool);
catch
    validChildren = allChildren;
end
set(figHandle,'Units','Pixels');
% set(validChildren,'Units','Pixels');

% types = get(allChildren,'Type');
% types = confirm_cell(types);
% uiMenuIndxsBool = cregexp(types,'uimenu');
% 
% %add all children except those of type "uimenu"
% validChildren = allChildren(~uiMenuIndxsBool);
% set(figHandle,'Units','Pixels');
% set(validChildren,'Units','Pixels');

% get the handles of the standalone figure
handles = guidata(figHandle);
if isempty(handles)
    handles = 'noguidata';
end

% get name of the standalone figure
figName = get(figHandle,'Name');
if isempty(figName) || strcmp(figName,' ')
    figName = ['tab ' num2str(tabNum)];
end

% create a tab
thisTabH = uitab(tabGroupH, ...
    'Title', figName, ...
    'UserData',tabNum, ...
    'Tag',get(figHandle,'Tag'), ... %make the original tabbedFig's tag this tab's tag
    'DeleteFcn',get(figHandle,'DeleteFcn'));%make the original tabbedFig's DeleteFcn this tab's DeleteFcn

% collect handles and tab info to tabbed gui's appdata 
guiAndTabInfo = getappdata(gcf,'guiAndTabInfo');
guiAndTabInfo(tabNum).handles = handles;
guiAndTabInfo(tabNum).tabHandles = thisTabH;
%remember the size of the original GUI and resize the tabbed GUI to this
%when tab is switched
guiAndTabInfo(tabNum).position =  get(figHandle,'Position');

%store info
setappdata(tabbedFig,'guiAndTabInfo', guiAndTabInfo);

% move objects from standalone figure to tab

% THERE APPEARS TO BE SOME PROBLEM PLOTTING LEGENDS (PARENT ASSOCIATION)
% Instead we get the vector of 'axes' types, as the highest level parents.
% set(validChildren(end),'Parent',thisTabH);
% set(validChildren,'Parent',thisTabH);
axesSubSet = findobj(validChildren,'Type','Axes');
set(axesSubSet,'Parent',thisTabH);                  

% close standalone figure since it has been "gutted" and placed onto a tab
delete(figHandle);
end

function update_guidata_and_resize(varargin)

if length(varargin) < 2
    return
end

% tabGroupH = varargin{1};
event_data = varargin{2};
if strcmp(event_data.EventName,'SelectionChange')
    curTabNum = get(event_data.NewValue,'UserData');
else
    return
end

guiAndTabInfo = getappdata(gcf,'guiAndTabInfo');
if isempty(guiAndTabInfo)
    return
end

%get handles of the children in the current tab
% handles = guiAndTabInfo(curTabNum).handles;
% 
% %update gui data with the handles of the children in the current tab
% if ~isempty(handles)
%     guidata(gcf, handles);
% end

newPos = calcNewPos(curTabNum, guiAndTabInfo);
set(gcf,'Position',newPos)

%force redraw
pause(0.01)
drawnow
end

function newPos = calcNewPos(curTabNum, guiAndTabInfo)
%update the size of the window to match the contents of the tab
%get position of gui when it opened as a standalone
figOrigPos = guiAndTabInfo(curTabNum).position;
newWidth = figOrigPos(3);
newHeight = 30 + figOrigPos(4); %assume tab is 30px tall

%ensure common units
set(gcf,'Units','Pixels');

%get current position
curFigPos = get(gcf,'Position');
curBottom = curFigPos(2);
curHeight = curFigPos(4);

% calculate new size
newBottom = curBottom + (curHeight-newHeight); %keep top left in place
newPos = [curFigPos(1), newBottom, newWidth, newHeight ];
end

function outCell = confirm_cell(inArg)
if ~iscell(inArg)
    outCell = {inArg};
else
    outCell = inArg;
end
end

function bool=cregexp(cellStrArray,pat)
%returns boolean array true at indices where pat is found in cellStrArray
cellStrArray = confirm_cell(cellStrArray);
bool = ~cellfun(@isempty,regexp(cellStrArray,pat));
end

function CenterWindow(hForeground, hBackground)
% centers gui with hForeground over gui hBackground
% hBackground is optional. If it's not included, hForeground is
% centered on the screen.
if nargin ==1
    %Center GUI Window
	origUnits = get(hForeground,'Units');
    set(hForeground,'Units','pixels');
    
    %get display size
    screenSize = get(0, 'ScreenSize');
    
    %calculate the center of the display
    newPos = get(hForeground, 'Position');
    newPos(1) = (screenSize(3)-newPos(3))/2;
    newPos(2) = (screenSize(4)-newPos(4))/2;
    
    %set new position of window
    set(hForeground, 'Position', newPos );
    set(hForeground,'Units',origUnits);
elseif nargin == 2
    % center hForeground over hBackground
	origUnitsF = get(hForeground,'Units');
	origUnitsB = get(hBackground,'Units');
    set(hForeground,'Units','pixels');
    set(hBackground,'Units','pixels');
    
    parentPos = get(hBackground, 'Position');
    
    %calculate the center of the parent, then offset by half the size of
    %hObject
    %     [left, bottom, width, height]
    parentCenter = [parentPos(1)+ parentPos(3)/2, parentPos(2)+ parentPos(4)/2];
    curPos = get(hForeground, 'Position');
    newPos(1) = parentCenter(1)- curPos(3)/2;
    newPos(2) = parentCenter(2)- curPos(4)/2;
    newPos(3) = curPos(3);
    newPos(4) = curPos(4);

    %set new position of foreground window
    set(hForeground, 'Position', newPos );
    
    set(hForeground,'Units',origUnitsF);
    set(hBackground,'Units',origUnitsB);
end
end