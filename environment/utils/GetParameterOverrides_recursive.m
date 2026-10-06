%% Recursive parameter search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [config] = GetParameterOverrides_recursive(defaultConfig,inputParameters)
% This function is designed to recursively search through a configuration
% structure/class for properties provided by a vector of parameter:value
% pairs.

% INPUTS:
% defaultConfig   - (cell/struct) defining original configuration.
% inputParameters - (cell array) of value pairs. 
% OUTPUTS:
% config    - The parameterised structure with 


% Default output state
config = defaultConfig;

% Sanity check #1 - Nothing to compare
if nargin < 2 
    return    
end

% If input is a structure, parse it
if isstruct(inputParameters)
    [inputParameters] = ParseStuctureToCellArray(inputParameters);
end

% Sanity check #2 - Check nested cell arrays
assert(iscell(inputParameters),'Entries must be provided as a cell array of parameter:value pairs');

% Sanity check #3 - Unloaded successfully
t = 1;
while numel(inputParameters) == 1 && t < 10
    inputParameters = inputParameters{:};
    t = t + 1;
end
assert(t < 10,'Failed to unload parameters.');

% Sanity check #4 - No inputs
if isempty(inputParameters) 
    return
end
    
% Sanity check #5 - Valid value/pair list
assert(mod(numel(inputParameters),2) == 0,'Un-even number of parameter:value pairs');

% ////////// BEGIN PARSING THE VECTOR AGAINST THE DATA STRUCTURE //////////
[config] = ScanDataStructure(config,inputParameters);
end

%% Convert structure to cell array of pairs
function [pairArray] = ParseStuctureToCellArray(inputStruct)

% Get the fieldnames
fn = fieldnames(inputStruct);
pairArray = cell(1,numel(fn)*2);

i = 1; j = 1;
while i ~= (numel(fn) + 1)
    pairArray{j} = fn{i};
    pairArray{j+1} = inputStruct.(fn{i});
    j = j + 2;
    i = i + 1;
end
end
%% Operations on structures 
function [config] = ScanDataStructure(config,inputParameters)
% Move through the properties of "config"
 
% Get the field names of the structure/class
fieldLabels = fieldnames(config); 

% Parse the top level structure
config = parseCellArray(config,inputParameters);

% Scan lower tier structures
for i = 1:length(fieldLabels)
    % Check if it's a structure or class
    if ~isstruct(config.(fieldLabels{i}))
        continue
    end
    % Conflict between "isobject" and type "sym"
    if isobject(config.(fieldLabels{i})) && ~isa(config.(fieldLabels{i}),'sym')
        continue
    end
    % Recursion
    config.(fieldLabels{i}) = ScanDataStructure(config.(fieldLabels{i}),inputParameters);
end
end
%% Parse the cell array against the data structure
function [config] = parseCellArray(config,inputParameters)

% Check the field against the parameter vector
for parameterIndex = 1:numel(inputParameters)
    % For each user-defined parameter
    givenParameter = inputParameters{parameterIndex};
    if ~ischar(givenParameter)
        continue
    end
    % If this is a structure
    if isstruct(config) && isfield(config,givenParameter)
        config.(givenParameter) = inputParameters{parameterIndex + 1};  % Make a substitution
    end
    % If this a class
    if ~isstruct(config) && isprop(config,givenParameter)
        config.(givenParameter) = inputParameters{parameterIndex + 1};   % Make a substitution
    end
end
end
