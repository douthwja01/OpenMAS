% CONVERT BODY AXIS RATES TO EULER RATES
function [eulerRates] = GetEularRatesFromOmega(eulerPosition,bodyAxisRates)
% This function calculates the equivalent eular rotational
% rates from defined body axis rates and a eular position.

% INPUT HANDLING
assert(numel(bodyAxisRates) == 3,'Incorrect number of body axis rates');
assert(numel(eulerPosition) == 3,'Incorrect number of Euler rotations');

% DEFINE THE CONVERSION MATRIX
conversionMatrix = eye(3);
conversionMatrix(1,2) = sin(eulerPosition(1))*tan(eulerPosition(2));
conversionMatrix(1,3) = cos(eulerPosition(1))*tan(eulerPosition(2));
conversionMatrix(2,2) = cos(eulerPosition(1));
conversionMatrix(2,3) = -sin(eulerPosition(1));
conversionMatrix(3,2) = sin(eulerPosition(1))*sec(eulerPosition(2));
conversionMatrix(3,3) = cos(eulerPosition(1))*sec(eulerPosition(2)); % VALIDATED

% DEFINE THE EULAR RATES FROM THE BODY AXIS RATES
eulerRates = conversionMatrix*bodyAxisRates;
end