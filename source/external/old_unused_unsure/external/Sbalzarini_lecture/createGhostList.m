% Code for Exercise 6 - Creation of ghost lists 

function ghostList = createGhostList(particles,numPerDim,cellList,numCells,cutoff)

dim = size(numCells,2);
numParticles=size(particles,1);
tempList=zeros(numParticles,1);
ghostList=[];
if dim==2
   
    for i=1:numParticles
        [I,J]=ind2sub(numCells,particles(i,dim+1));
        if (I==1 | I==numCells(1)| J==1 | J==numCells(2));
            tempList(i)=1;
            
        end

    end
    
    ghostList=find(tempList);
    
elseif dim==3
    
    for i=1:numParticles
        [I,J,K]=ind2sub(numCells,particles(i,dim+1));
        if (I==1 | I==numCells(1)| J==1 | J==numCells(2) | K==1 | K==numCells(3));
            tempList(i)=1;
        end

    end
    
    ghostList=find(tempList);
    
end
