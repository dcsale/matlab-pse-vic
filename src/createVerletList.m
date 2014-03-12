function verletList = createVerletList(particleMat,cellList,numCells,cutoff)
%%%%%%%%%%%%%%%%%%%%%%%% 
% Creation of Verlet lists without symmetry
%%%%%%%%%%%%%%%%%%%%%%%% 
% Input
% particleMat:  (numParticles x dim+1)-Matrix of particle positions and cell index
% cellList:     Matlab cell structure of length prod(numCells) (number of cells)
% numCells:     (dim x 1)- Vector that contains the number of cells per
%               dimension
% cutoff:       Scalar distance cutoff that defines the neighborhood. It should be
%               equivalent to the side length of a cell
% 
% Output
% verletList:   Matlab cell structure of length numParticles. Each element
%               contains the indices of the particles, i.e. the row numbers
%               in the particleMat matrix
%
% function verletList = createVerletList(particleMat,cellList,numCells,cutoff)

[numParticles,dim] = size(particleMat);
dim = dim-1;

% Initialize Verlet list
verletList = cell(numParticles,1);

% Iterate over all particles and build Verlet list
for n=1:numParticles
    % collect all surrounding cell indices
    currCellInd = particleMat(n,dim+1);
    adjCellIndices = adjacentCells(currCellInd,numCells);
    allCellIndices = [currCellInd,adjCellIndices]';
    
    % collect all particle ID's
    currPartIndices=cell2mat(cellList(allCellIndices));
    
    currVec=[];
    numRows = length(currPartIndices);
    
    % iterate over all collected particles around the current particle
    tempMat = repmat(particleMat(n,1:dim),numRows,1);
    temp = currPartIndices(find(sum(((tempMat-particleMat(currPartIndices,1:dim)).^2),2)<=cutoff^2));
    verletList{n} = temp(find(temp~=n));
end

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

% Plot the selected particle and its neighbors

figure
plot(particleMat(randParticle,1),particleMat(randParticle,2),'m.','MarkerSize',10)
xlim([lBounds uBounds])
ylim([lBounds uBounds])
zlim([lBounds uBounds])
set(gca,'XTick',[lBounds(1):cutoff:uBounds(1)])
set(gca,'YTick',[lBounds(1):cutoff:uBounds(1)])
set(gca,'ZTick',[lBounds(1):cutoff:uBounds(1)])
grid on
box on
xlabel('x')
ylabel('y')
zlabel('z')
hold on

verletParts = verletList{randParticle};

% Plot neighbors
plot(particleMat(verletParts,1),particleMat(verletParts,2),'r.','MarkerSize',20)

% Plot transparent cutoff sphere
[x, y, z] = ellipsoid(particleMat(randParticle,1),particleMat(randParticle,2),0,cutoff,cutoff,cutoff,30);
surf(x, y, z)
shading interp
colormap(gray);
alpha(.4)



