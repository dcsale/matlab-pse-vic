function varargout = global_inds(d, varargin)
%GLOBAL_INDS Returns the global indices for all
%processors in the map of distributed array D.
%   GLOBAL_INDS(D, DIM) Returns global indices  of the 
%       distributed array D for all processors in the 
%       specified dimension, DIM.
%   GLOBAL_INDS(D) Returns global indices of the 
%       distributed array D for all processors in 
%       all dimensions of D.
%   For each dimension, the following is the format of the indices
%   returned:
%       Matrix of size NUM_PROCS_IN_GRIDxMAX_LOCAL_INDS. Each line of the returned
%       matrix M, M(i,:) contains the follwing information 
%       [PROCESSOR_RANK IND1 IND2 ... INDn]
%       To ensure that all rows in the return index are the same, the
%       indices matrix is appended with extra zeros where there are not
%       enough indices.
%
% Author:   Nadya Travinin

%dimension of the distributed object
dim = length(d.pitfalls);

%processor grid on which the object is distributed
grid = d.map.grid;

len = size(grid);
if length(len)<dim
    for i = (length(len)+1):dim
        len(i) = 1;
    end
end

if dim==2
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            local_falls = get_local_falls(d.pitfalls, grid, grid(i,j));
            global_ind{i,j} = get_global_ind(local_falls, len);
        end
    end
elseif dim==3
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                local_falls = get_local_falls(d.pitfalls, grid, grid(i,j,k));
                global_ind{i,j,k} = get_global_ind(local_falls, len);
            end
        end
    end
elseif dim==4
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                for m = 1:len(4)
                    local_falls = get_local_falls(d.pitfalls, grid, grid(i,j,k,m));
                    global_ind{i,j,k,m} = get_global_ind(local_falls, len);
                end
            end
        end
    end
else
    error('GLOBAL_BLOCK_RANGES: Only objects up to 4-D are supported');
end



%----------GLOBAL_BLOCK_RANGE--------------------------------------
if length(varargin)==1
    dims = varargin{1}; %desired index dimension is passed in
else
    dims = 1: length(d.size); %generate ranges for all dimensions
end

my_inds = d.global_ind;
s = d.size;
for i = 1:length(dims)
    clear temp;
    num_procs = prod(len); %total number of grid processors
    proc_count = 1; %keep track of the number of processors
    
    if dim==2 %2D array
        for g1 = 1:len(1) %grid cols
            for g2=1:len(2) %grid rows
                curr_inds = global_ind{g1,g2}; 
                if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                    ind_len = 1+s(dims(i));
                    temp(proc_count,1:ind_len) = [grid(g1,g2) 1:s(dims(i))]; %block range is the entire dimension
                else %block range is not the entire range
                    dim_inds = curr_inds{dims(i)};
                    ind_len = 1+length(dim_inds);
                    temp(proc_count,1:ind_len) = [grid(g1,g2) dim_inds];
                end
                proc_count = proc_count+1;
            end
        end
    elseif dim==3 %3D array
        for g1 = 1:len(1)
            for g2=1:len(2)
                for g3 = 1:len(3)
                    curr_inds = global_ind{g1,g2,g3};
                    if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                        ind_len = 1+s(dims(i));
                        temp(proc_count,1:ind_len) = [grid(g1,g2,g3) 1:s(dims(i))]; %block range is the entire dimension
                    else %block range is not the entire range
                        dim_inds = curr_inds{dims(i)};
                        ind_len = 1+length(dim_inds);
                        temp(proc_count,1:ind_len) = [grid(g1,g2,g3) dim_inds];
                    end
                    proc_count = proc_count+1;
                end
            end
        end
    elseif dim==4 %4D array
        for g1 = 1:len(1)
            for g2=1:len(2)
                for g3 = 1:len(3)
                    for g4 = 1:len(4)
                        curr_inds = global_ind{g1,g2,g3,g4};
                        if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                            ind_len = 1+s(dims(i));
                            temp(proc_count,1:ind_len) = [grid(g1,g2,g3,g4) 1:s(dims(i))]; %block range is the entire dimension
                        else %block range is not the entire range
                            dim_inds = curr_inds{dims(i)};
                            ind_len = 1+length(dim_inds);
                            temp(proc_count,1:ind_len) = [grid(g1,g2,g3,g4) dim_inds];
                        end
                        proc_count = proc_count+1;
                    end
                end
            end
        end         
    end %dim
    varargout{i} = temp;
end
%------------------------------------------------------------------

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
