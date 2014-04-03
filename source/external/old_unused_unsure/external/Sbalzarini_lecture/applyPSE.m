%%%%%%%%%%%%%%%%%%%%%%%%
% 2D PSE Operator
%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% particleMat:  (numParticles x (dim+1+numStren))-Matrix of particle positions, cell
%               indices and particle strengths 
% verletList:   Verlet list of particles
% epsilon:      Kernel parameter epsilon (standard deviation)
% numStren:     Number of strengths a particle carries
%
% Output
% pseSum:       ((numParticles x numStren)-Matrix of updated particle strengths
%
% function pseSum = applyPSE(particleMat,verletList,epsilon,numStren)

function pseSum = applyPSE(particleMat,verletList,epsilon,numStren)

dim = size(particleMat,2)-1-numStren;
numParticles = length(verletList);
pseSum = zeros(numParticles,numStren);

for i=1:numParticles

    % collect neighboring particles from the Verlet list
    neighParticleMat = particleMat(verletList{i},:);
    pseSum(i,1:2) = 0;

    if ~isempty(neighParticleMat)
        % compute distances
        neighVecs=neighParticleMat(:,1:dim)-repmat(particleMat(i,1:dim),size(neighParticleMat,1),1);
        particleDists=sqrt(sum(neighVecs.*neighVecs,2));

        % Compute summation of the PSE operator
        for j=1:numStren
            pseSum(i,j) = sum((neighParticleMat(:,dim+1+j)-particleMat(i,dim+1+j)).*etaKernel(particleDists,epsilon,dim));
        end
    else
        pseSum(i,1) = 0;
    end


end

