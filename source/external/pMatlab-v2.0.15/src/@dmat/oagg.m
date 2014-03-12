function mat = oagg(d)
%AGG Aggregates the parts of a distributed matrix on the leader processor.
%    AGG(D) aggregates the parts of a distributed matrix 
%    into a whole and returns it as a regular matrix.
%    If the current processor is the LEADER, returns 
%    the aggreagated matrix, otherwise, returns the 
%    local part.
%    This functions increments GLOBAL message TAG.
%   
%    NOTE: Currently, it doesn't matter whether or not the leader is in the
%    map - the global matrix is returned on the leader regardless. 
%
% Author:   Nadya Travinin
           
%globals
global pMATLAB;

%increment tag
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

if pMATLAB.my_rank == pMATLAB.leader
    if d.dim==2
        dim = size(d.map.grid);
        %dim(1) - number of grid rows, dim(2) - number of grid cols
        for i = 1:dim(1)
            for j = 1:dim(2)
                if (pMATLAB.my_rank==d.map.grid(i,j))
                    temp_mat{i,j} = d.local;
                else
                    temp_mat{i,j} = MPI_Recv(d.map.grid(i,j), pMATLAB.tag, pMATLAB.comm);
                end
            end
        end
        
        %reconstruct the matrix from the local pieces
        %this is a NO-OP for block distributions since the data does not
        mat = reconstruct(d.pitfalls,  d.map.grid, temp_mat, d.size); 
    elseif d.dim==3
        dim = size(d.map.grid);
        if length(dim)<d.dim
            for i = (length(dim)+1):d.dim
                dim(i) = 1;
            end
        end
        for i = 1:dim(1)
            for j = 1:dim(2)
                for k = 1:dim(3)
                    if (pMATLAB.my_rank==d.map.grid(i,j, k))
                        temp_mat{i,j,k} = d.local;
                    else
                        temp_mat{i,j,k} = MPI_Recv(d.map.grid(i,j,k), pMATLAB.tag, pMATLAB.comm);
                    end
                end
            end
        end
    
        %reconstruct the matrix from the local pieces
        %this is a NO-OP for block distributions since the data does not
        mat = reconstruct(d.pitfalls,  d.map.grid, temp_mat, d.size); 
        
    elseif d.dim==4
        dim = size(d.map.grid);
        if length(dim)<d.dim
            for i = (length(dim)+1):d.dim
                dim(i) = 1;
            end
        end
        for i = 1:dim(1)
            for j = 1:dim(2)
                for k = 1:dim(3)
                    for m = 1:dim(4)
                        if (pMATLAB.my_rank==d.map.grid(i,j,k,m))
                            temp_mat{i,j,k,m} = d.local;
                        else
                            temp_mat{i,j,k,m} = MPI_Recv(d.map.grid(i,j,k,m), pMATLAB.tag, pMATLAB.comm);
                        end
                    end
                end
            end
        end
    
        %reconstruct the matrix from the local pieces
        %this is a NO-OP for block distributions since the data does not
        mat = reconstruct(d.pitfalls,  d.map.grid, temp_mat, d.size); 
    else
        error('DMAT/AGG: Only up to 4-D objects currently supported');
    end
else %my_rank ~= leader
    %send local data to the leader regardless of the matrix dimension
    if inmap(d.map, pMATLAB.my_rank) %only send data if processor is in the map
        MPI_Send(pMATLAB.leader, pMATLAB.tag, pMATLAB.comm, d.local);
    end
    mat = d.local;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMatlab: Parallel Matlab Toolbox
% Software Engineer: Ms. Nadya Travinin (nt@ll.mit.edu)
% Architect:      Dr. Jeremy Kepner (kepner@ll.mit.edu)
% MIT Lincoln Laboratory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2005, Massachusetts Institute of Technology All rights 
% reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are 
% met:
%      * Redistributions of source code must retain the above copyright 
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright 
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the
%        distribution.
%      * Neither the name of the Massachusetts Institute of Technology nor 
%        the names of its contributors may be used to endorse or promote 
%        products derived from this software without specific prior written 
%        permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
