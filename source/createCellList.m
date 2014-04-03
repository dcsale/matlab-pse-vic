function [particleMat, cellList, numCells] = createCellList(xp, MESH)
%%%%%%%%%%%%%%%%%%%%%%%% 
% Creation of cell lists 
%%%%%%%%%%%%%%%%%%%%%%%% 
% Input
% xp:           (numParticles x dim) Matrix of particle positions 
% xp:           (numParticles x dim) Matrix of particle positions 
% lBounds:      Scalar lower bound on all particle positions  
% uBounds:      Scalar upper bound on all particle positions
% 
% Output
% particleMat:  (numParticles x dim+1)-Matrix that contains the particle
%               positions and the cell index it belongs to
% cellList:     Matlab cell structure of length number of overall cells (prod(numCells)) 
% numCells:     (dim x 1) Vector that contains the number of cells per dimension 

xp = xp';


% cellSide:     Scalar (or Vector???) value of the cell's side length 
% cellSide = MESH.dx(1);
[numParticles, dim] = size(xp);
switch dim
    case 1
        cellSide = MESH.dx(1);
    case 2
        cellSide = [MESH.dx(1), MESH.dx(2)];
    case 3
        cellSide = [MESH.dx(1), MESH.dx(2), MESH.dx(3)];
    otherwise
        error('error [createCellList.m]: dim > 3 not recognized')
end
 
% domainSize = (uBounds - MESH.xmin)*ones(dim,1);
domainSize = MESH.xmax(1:dim) - MESH.xmin(1:dim);

% determine how many cells are needed per dimension
% numCells = max(1,floor(domainSize./cellSide)); % original for scalar cellSide
numCells = max(1, floor(domainSize./cellSide)); % modified for vector cellSide
indParticlePos = zeros(numParticles, dim);

parfor i=1:dim
    % shift to zero and discretize
    indParticlePos(:,i) = floor( (xp(:,i) - min(xp(:,i))) ./ (domainSize(i)-min(xp(:,i))).*numCells(i) ) + 1;
end

% add cell index in the particle representation
indices = zeros(numParticles, 1);
switch dim
    case 1
        indices = indParticlePos;
        
    case 2
        ix = indParticlePos(:,1);
        iy = indParticlePos(:,2);
        parfor i = 1:numParticles
            % Use the convention of sub2ind to map the dim x1 indices to a single integer
%             indices(i) = sub2ind(numCells, indParticlePos(i,1), indParticlePos(i,2));
            indices(i) = sub2ind(numCells, ix(i), iy(i));
        end
        
    case 3
        ix = indParticlePos(:,1);
        iy = indParticlePos(:,2);
        iz = indParticlePos(:,3);
        parfor i = 1:numParticles
            % Use the convention of sub2ind to map the dim x1 indices to a single integer
%             indices(i) = sub2ind(numCells, indParticlePos(i,1), indParticlePos(i,2), indParticlePos(i,3));
            indices(i) = sub2ind(numCells, ix(i), iy(i), iz(i));
        end
        
    otherwise
        error('dim > 3 not supported')
end

% Create a Matlab cell structure
particleMat = [xp, indices];
maxCellNum  = max( particleMat(:,dim+1) );
cellList    = cell(prod(numCells), 1);

cellInd = particleMat(:,dim+1);
parfor c = 1:maxCellNum
%     currInd     = find(particleMat(:,dim+1)==c);
    currInd     = find(cellInd==c);
    cellList{c} = currInd;    
end

% plot_random = true;
plot_random = false;
if plot_random
    test_CellList(particleMat, cellList, numCells, dim, MESH);
end

end % function createCellList


function test_CellList(particleMat, cellList, numCells, dim, MESH)
%% Test and plot the implementation

% Determine the particles in a random cell
maxCellInd = max(particleMat(:,dim+1));
randCell   = ceil(rand*maxCellInd);

% show the randomly chosen 2D/3D cell index
% [I,J] = ind2sub(numCells,randCell);
[I, J, K] = ind2sub(numCells,randCell);


% There are two possibilities to extract the particles contained in the
% randomly chosen cell
% (i) Collect the particles from the particle matrix
% randCellIndices=find(particleMat(:,dim+1)==randCell);
% (ii) Collect the particles from the cell list
randCellIndices = cellList{randCell};

figure
switch dim
    case 2
        title(['Displaying the cell with indices I =',int2str(I),', J =',int2str(J)])
        plot(particleMat(randCellIndices,1),particleMat(randCellIndices,2),'c.')
        xlim([MESH.xmin(1) MESH.xmax(1)])
        ylim([MESH.xmin(2) MESH.xmax(2)])
%         set(gca, 'XTick', MESH.xmin(1):cutoff:MESH.xmax(1))
%         set(gca, 'YTick', MESH.xmin(2):cutoff:MESH.xmax(2))
        set(gca, 'XTick', MESH.xmin(1):MESH.dx(1):MESH.xmax(1))
        set(gca, 'YTick', MESH.xmin(2):MESH.dx(2):MESH.xmax(2))
        xlabel('x')
        ylabel('y')

    case 3
        title(['Displaying the cell with indices I =' int2str(I) ', J =',int2str(J) ', K =' int2str(K)])
        plot(particleMat(randCellIndices,1),particleMat(randCellIndices,2),'c.')
        
end
grid on
box on

hold on
switch dim
    case 2
        plot(particleMat(randCellIndices,1), particleMat(randCellIndices,2), 'bo')
        
    case 3
%         plot(particleMat(randCellIndices,1),particleMat(randCellIndices,2),'bo')
end


%%%%%%%%%%%%%%%%%%%
% Check adjacency implementation by determining the adjacent cells of the
% randomly chosen cell

adjCellIndices = adjacentCells(randCell, numCells)';
adjPartIndices = cell2mat(cellList(adjCellIndices));

plot(particleMat(adjPartIndices,1), particleMat(adjPartIndices,2), 'y.')

% %%%%%%%%%%%%%%%%%%%%%
% % Create Verlet list
% verletList = createVerletList(particleMat,cellList,numCells,cutoff);
% 
% % Choose the first particle in the previously chosen random cell
% randParticle = cellList{randCell};
% randParticle = randParticle(1);

end % function test_CellList



