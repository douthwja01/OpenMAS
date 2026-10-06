% PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
function [config] = GetParameterOverrides(defaultConfig,inputParameters)
% This function is designed to compare a vector of {'input',parameter} pairs
% against a structure(or object) with labelled fields.

% INPUTS:
% defaultConfig     - The reference structure/class with available fields.
% inputParameters   - The cell array of string, parameter pairs or
%                     comparative structure.
% OUTPUT:
% config            - The structure/class with the fields overitten. 

% Input sanity check #1
if nargin < 2 || numel(inputParameters) == 0
    config = defaultConfig;
    return
end

% Input sanity check #2 - Check nested cell arrays
if iscell(inputParameters)
    t = 1;
    while numel(inputParameters) == 1 && t < 10
        inputParameters = inputParameters{:};
        t = t + 1;
    end
    assert(t ~= 10,'Failed to unload parameters.');
end

% Input-dependant behaviour
switch true
    case iscell(inputParameters)
        assert(mod(numel(inputParameters),2) == 0,'Entries are not of list of parameter:value pairs');
        config = parseCellArray(defaultConfig,inputParameters);
        return
    case isstruct(inputParameters)
        config = parseStructure(defaultConfig,inputParameters);
        return
    otherwise
        error('Please provide a list of parameter:value pairs or a structure.');
end
end

% Parse the input parameters against a list of parameter:value cell array
function [config] = parseCellArray(config,inputParameters)

% Move through the parameter pairs (i.e. ,'radius',1,...)
for parameterIndex = 1:numel(inputParameters)
    % For each user-defined parameter
    givenParameter = inputParameters{parameterIndex};
    if ~ischar(givenParameter)
        continue
    end
    
    % IF THE OBJECT HAS A PROPERTY BY THAT NAME
    if isstruct(config) && isfield(config,givenParameter)
        config.(givenParameter) = inputParameters{parameterIndex + 1};  % Make a substitution
    end
    if ~isstruct(config) && isprop(config,givenParameter)
        config.(givenParameter) = inputParameters{parameterIndex + 1};   % Make a substitution
    end
end
end
% Parse the input parameters against a structure
function [config] = parseStructure(config,inputStructure)

% Move through the fields of the second structure and override associated
% fields.
fields = fieldnames(inputStructure);
for f = 1:numel(fields)
    % Config is a struct
    if isstruct(config) && isfield(config,fields{f})
        config.(fields{f}) = inputStructure.(fields{f});        % Make a substitution
        continue
    end
    % Config is a class
    if ~isstruct(config) && isprop(config,fields{f})
        config.(fields{f}) = inputStructure.(fields{f});        % Make a substitution
        continue
    end
end
end