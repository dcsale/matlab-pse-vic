%%%%%%%%%%%%%%%%%%%%%%%% 
% Find adjacent cells given a cell index and the
% size of the cells
%%%%%%%%%%%%%%%%%%%%%%%% 
% Input
% cellIndex:       Integer representing cell index
% numCells:        (dim x 1)-Vector with number of cells 
%
% 
% Output
% adjCellIndices:  Vector of adjacent cell indices 
%
% function adjCellIndices = adjacentCells(cellIndex,numCells)

function adjCellIndices = adjacentCells(cellIndex,numCells)

dim=size(numCells,1);
adjCellIndices=[];

if dim==3
    
    [I,J,K]=ind2sub(numCells,cellIndex);
    
    iIndices=[max(I-1,1):1:min(I+1,numCells(1))];
    jIndices=[max(J-1,1):1:min(J+1,numCells(2))];
    kIndices=[max(K-1,1):1:min(K+1,numCells(3))];
    
    for i=iIndices
        for j=jIndices
            for k=kIndices
                if any([i,j,k]-[I,J,K])
                    adjCellIndices=[adjCellIndices,sub2ind(numCells,i,j,k)];
                end
            end
        end
    end
        
elseif dim==2
    
    [I,J]=ind2sub(numCells,cellIndex);
    
    iIndices=[max(I-1,1):1:min(I+1,numCells(1))];
    jIndices=[max(J-1,1):1:min(J+1,numCells(2))];
        
    for i=iIndices
        for j=jIndices
            if any([i,j]-[I,J])
                    adjCellIndices=[adjCellIndices,sub2ind(numCells,i,j)];
            end
        end
    end
end