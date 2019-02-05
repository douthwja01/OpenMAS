% CHECK THE AGENT VELOCITY IS INSIDE THE VO
function [flag] = isInsideVO(point,VO,isSymmetric)
% determine if the point p is inside the given VO
% angle<halfOpenAngle; here numericalTolerance is an error tolarance

if nargin < 3
    isSymmetric = 1;
end

flag = 0;
VOtolerance = 1E-8;
% Vector to point
candidateVector = point - VO.apex;

% Evaluate based on symmetry of the VO
if isSymmetric
    
    VOprojection = norm(candidateVector)*cos(VO.openAngle/2);
    
    candProjection = VO.axisUnit'*candidateVector;
    projDiff = (candProjection - VOprojection);
    if projDiff > VOtolerance
        flag = 1;
    end
else
    %                 VOprojection = norm(candidateVector)*cos(VO.openAngle/2);
end

end