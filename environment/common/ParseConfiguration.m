% PARSE A GENERIC INPUT SET AGAINST A DEFAULT CONFIG STRUCTURE
function [config] = ParseConfiguration(defaultConfig,inputParameters)
% This function is designed to compare a vector of {'input',parameter} pairs
% against a structure(or object) with labelled fields.

% INPUTS:
% defaultConfig - The reference structure/class with available fields.
% inputParameters - The cell array of string, parameter pairs.
% OUTPUT:
% config - The structure/class with the fields overitten. 

% Input sanity check
assert(nargin == 2,'Parser must be provided with an object/structure and a cell parameter vector');
assert(iscell(inputParameters),'The input vector must be a cell array of parameter; value pairs.');
assert(mod(numel(inputParameters),2) == 0,' Please provide list of parameter:value pairs');

% MOVE THROUGH THE PARAMETER PAIRS (i.e. ,'radius',1,...)
for parameterIndex = 1:numel(inputParameters)
    % FOR EACH USER-DEFINED PARAMETER
    givenParameter = inputParameters{parameterIndex};
    if ~ischar(givenParameter)
        continue
    else
        % IF THE OBJECT HAS A PROPERTY BY THAT NAME
        if isstruct(defaultConfig) && isfield(defaultConfig,givenParameter)
            defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};  % Make a substitution
        end
        if ~isstruct(defaultConfig) && isprop(defaultConfig,givenParameter)
            defaultConfig.(givenParameter) = inputParameters{parameterIndex + 1};   % Make a substitution
        end
    end
end
config = defaultConfig;
end