function [i,j] = find(x)
%FIND Find indices of nonzero elements.
%   [I,J] = FIND(X) returns the row and column indices of nonzero elements
%   of the distributed matrix X. 
%   NOTE: Currently supports only [i,j] = find(x) calling convention.
%         Only works on 2D arrays.
%
% Author:   Nadya Travinin

%global variables
global pMATLAB;


%increment tag
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));    

if inmap(x.map, pMATLAB.my_rank)
    [local_i, local_j] = find(x.local);
    if isa(x.global_ind{1}, 'char')
        if x.global_ind{1} == ':'
            x.global_ind{1} = [1:x.size(1)];
        end
    end
    if isa(x.global_ind{2}, 'char')
        if x.global_ind{2} == ':'
            x.global_ind{2} = [1:x.size(2)];
        end
    end

    % When a processor is allocated a single row (column), x.global_ind{1}
    % (x.global_ind{2}) contain a 1x1 matrix.  Since local_i and local_j
    % are
    % column vectors, global_i (global_j) will contain a column vector.
    % Transpose local_i and local_j so that global_i and global_j are row
    % vectors.
    local_i = local_i';
    local_j = local_j';
    global_i = x.global_ind{1}(local_i);
    global_j = x.global_ind{2}(local_j);
    
    data.i = global_i;
    data.j = global_j;

    grid_size = size(x.map.grid);
    %grid_size(1) - number of grid rows, grid_size(2) - number of grid cols
    %send local finds to everyone
    for d1 = 1:grid_size(1)
        for d2 = 1:grid_size(2)
            if (pMATLAB.my_rank~=x.map.grid(d1, d2))
                MPI_Send(x.map.grid(d1, d2), pMATLAB.tag, pMATLAB.comm, data);
            end
        end
    end
    
    %receive finds from everyone
    for d1 = 1:grid_size(1)
        for d2 = 1:grid_size(2)
            if (pMATLAB.my_rank~=x.map.grid(d1, d2))
                temp{d1, d2} = MPI_Recv(x.map.grid(d1, d2), pMATLAB.tag, pMATLAB.comm);
            else
                temp{d1, d2} = data;
            end
        end
    end

    i = [];
    j = [];

    for d2 = 1:grid_size(2) %grid cols
        for d1 = 1:grid_size(1) %grid rows
            if (pMATLAB.my_rank~=x.map.grid(d1, d2))
                if ~isempty(temp{d1,d2}.i)
                    i = [i temp{d1,d2}.i];
                    j = [j temp{d1,d2}.j];
                end
            else
                if ~isempty(data.i)
                    i = [i data.i];
                    j = [j data.j];
                end
            end
        end
    end
    %transpose outputs since MATLAB find returns column vectors
    i = i';
    j = j';
else
    i = [];
    j = [];
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
