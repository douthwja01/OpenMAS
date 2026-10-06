% CONVERT EULER RATES TO BODY AXIS RATES
function [bodyAxisRates] = GetOmegaFromEulerRates(eulerPosition,eulerRates)
% This function calculates the conversion between the eular
% rates and the body axis rates.

% INPUT HANDLING
assert(numel(eulerRates) == 3,'Incorrect number of Euler rates');
assert(numel(eulerPosition) == 3,'Incorrect number of Euler rotations');

% DEFINE THE CONVERSION MATRIX
conversionMatrix = eye(3,3);
conversionMatrix(1,3) = -sin(eulerPosition(2));
conversionMatrix(2,2) =  cos(eulerPosition(1));
conversionMatrix(2,3) =  cos(eulerPosition(2))*sin(eulerPosition(1));
conversionMatrix(3,2) = -sin(eulerPosition(1));
conversionMatrix(3,3) =  cos(eulerPosition(2))*cos(eulerPosition(1)); % VALIDATED
% MAP THE EULAR ROTATION RATES TO THE BODY AXIS ROTATIONAL RATES
bodyAxisRates = conversionMatrix*eulerRates;
end