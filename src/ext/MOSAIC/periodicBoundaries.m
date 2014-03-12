%%%%%%%%%%%%%%%%%%%%%%%%
% Periodic boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% strengthMat:   (sqrt(numParticles) x sqrt(numParticles))-Matrix of particle 
%                strengths 
% h:             grid spacing
% cutoff:        Cutoff value for the neighborhood
% 
% Output
% strengthBound: ((numParticles x dim+2)-Matrix of particle positions, cell
%                indices and updated particle strengths
%
% function strengthBound = periodicBoundaries(strengthMat,h,cutoff)

function strengthBound = periodicBoundaries(strengthMat,h,cutoff)

bW = round(cutoff/h);
[rows,cols] = size(strengthMat);

strengthBound = strengthMat;

% permute left/right boundary
strengthBound(:,cols-(bW-1:-1:0)) = strengthMat(:,bW+1:2*bW);
strengthBound(:,1:bW) = strengthMat(:,cols-(2*bW-1:-1:bW));

% permute upper/lower boundary
strengthBound(rows-(bW-1:-1:0),:) = strengthMat(bW+1:2*bW,:);
strengthBound(1:bW,:) = strengthMat(rows-(2*bW-1:-1:bW),:);

% permute the remaining corners up/down
strengthBound(1:bW,1:bW) = strengthMat(rows-(2*bW-1:-1:bW),cols-(2*bW-1:-1:bW));
strengthBound(rows-(bW-1:-1:0),cols-(bW-1:-1:0)) = strengthMat(bW+1:2*bW,bW+1:2*bW);

strengthBound(1:bW,cols-(bW-1:-1:0)) = strengthMat(rows-(2*bW-1:-1:bW),bW+1:2*bW);
strengthBound(rows-(bW-1:-1:0),1:bW) = strengthMat(bW+1:2*bW,cols-(2*bW-1:-1:bW));

