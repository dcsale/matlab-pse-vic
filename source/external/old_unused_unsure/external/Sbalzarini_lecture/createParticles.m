%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of a uniformly random or grid-based
% distribution of particles 
%%%%%%%%%%%%%%%%%%%%%%%% 
% Input
% numParticles: Number of particles
% dim:          Dimension (dim==1 or dim==2)
% lBounds:      Scalar lower bound on all particle positions  
% uBounds:      Scalar upper bound on all particle positions
% kind:         Integer defining the kind of sampling (k==0 random, k==1 grid with particles 
%               in the center of a cell.)
%
% Output
% particlePos:    (numParticles x dim)-Matrix of particle positions   
%
% function particlePos =
% createParticles(numParticles,dim,lBounds,uBounds,kind)
function particlePos = createParticles(numParticles,dim,lBounds,uBounds,kind)

particlePos = [];

if kind == 0
    
    % Particles that are distributed uniformly at random 
    particlePos = lBounds + (uBounds-lBounds).* rand(numParticles,dim);
    
elseif kind == 1

    % Particles that are distributed on a grid 
    numPerDim=round(numParticles^(1/dim));
    cellWidth = (uBounds-lBounds)/(2*numPerDim);
    
    if dim==1
        
        particlePos = linspace(lBounds+cellWidth,uBounds-cellWidth,numPerDim);
        
    elseif dim==2
        
        xVec = linspace(lBounds+cellWidth,uBounds-cellWidth,numPerDim);
        yVec = linspace(lBounds+cellWidth,uBounds-cellWidth,numPerDim);
        [X,Y]=meshgrid(xVec,yVec);
        particlePos=[X(:),Y(:)]; 
        
    end
end
