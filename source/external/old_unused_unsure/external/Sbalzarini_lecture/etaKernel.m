%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel function for PSE operator
%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% rVec:         (numNeighbors x 1)-Vector of inter-particle distances
% epsilon:      Kernel parameter epsilon (standard deviation)
% dim:          dimension
%
% Output
% etaVals:      (numNeighbors x 1)-Vector of eta values
%
% function etaVals = etaKernel(r,epsilon,dim)

function etaVals = etaKernel(rVec,epsilon,dim)

switch dim
    case 1
          etaVals = 1/(2*epsilon*sqrt(pi)) * exp(-rVec.^2./(4*epsilon^2));
    case 2
          etaVals = 4/(epsilon^2*pi) * exp(-rVec.^2./(epsilon^2));
    case 3
          etaVals = ones(size(rVec))*(15/pi^2)./(abs(rVec).^10+ones(size(rVec)));    
end