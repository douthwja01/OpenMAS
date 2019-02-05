% CONVERT VECTOR TO SKEW-SYMMETRIC MATRIX
function [outputMatrix] = getSkewSymmetricMatrix(inputVector)
% This function generates a skew-symmetric for the computation
% of the vector cross-product.
% INPUTS:
% inputVector - The original 3D vector
% OUTPUT:
% outputMatrix - The equivalent skew-symmetric matrix

if length(inputVector) ~= 3
    warning('The input vector must be three dimensional.');
    return
end

outputMatrix = zeros(3,3);
outputMatrix(1,2) = -inputVector(3);
outputMatrix(1,3) =  inputVector(2);
outputMatrix(2,1) =  inputVector(3);
outputMatrix(2,3) = -inputVector(1);
outputMatrix(3,1) = -inputVector(2);
outputMatrix(3,2) =  inputVector(1); % Arrange the components
end