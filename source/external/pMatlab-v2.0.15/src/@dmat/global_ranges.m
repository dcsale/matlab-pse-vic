function varargout = global_ranges(d, varargin)
%GLOBAL_RANGES Returns the ranges of global indices for all
%   processors in the map of distributed array D. Returns the same range if D
%   is block distributed, returns subranges for (block) cyclic distirbutions.
%   GLOBAL_RANGES(D, DIM) Returns the global index ranges of the 
%       distributed array D for all processors in the 
%       specified dimension, DIM.
%   GLOBAL_RANGES(D) Returns the global index range of the 
%       distributed array D for all processors in 
%       all dimensions of D.
%   For each dimension, the following is the format of the indices
%   returned:
%       Matrix of size NUM_PROCS_IN_GRIDx(NUM BLOCK BOUNDARIES). Each line of
%       the returned matrix M, M(i,:) contains the follwing information 
%       [PROCESSOR_RANK START_INDEX_1 END_INDEX_1 START_INDEX_2 END_INDEX_2 ...]
%   NOTE:
%       If processors in the dimension have different number of blocks, the
%       block boundaries are padded with zeros for the processors that have
%       fewer blocks.
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
            local_falls{i,j} = get_local_falls(d.pitfalls, grid, grid(i,j));
            global_ind{i,j} = get_global_ind(local_falls{i,j}, len);
        end
    end
elseif dim==3
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                local_falls{i,j,k} = get_local_falls(d.pitfalls, grid, grid(i,j,k));
                global_ind{i,j,k} = get_global_ind(local_falls{i,j,k}, len);
            end
        end
    end
elseif dim==4
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                for m = 1:len(4)
                    local_falls{i,j,k,m} = get_local_falls(d.pitfalls, grid, grid(i,j,k,m));
                    global_ind{i,j,k,m} = get_global_ind(local_falls{i,j,k,m}, len);
                end
            end
        end
    end
else
    error('GLOBAL_RANGES: Only objects up to 4-D are supported');
end



%----------GLOBAL_RANGE--------------------------------------
if length(varargin)==1
    dims = varargin{1}; %desired index dimension is passed in
else
    dims = 1: length(d.size); %generate ranges for all dimensions
end

my_inds = d.global_ind;
my_falls = d.falls;
s = d.size;

for i = 1:length(dims)
    num_procs = prod(len); %total number of grid processors
    num_bounds = my_falls(dims(i)).n*2; %num bounds is the number of line segments*2 (beginning and end of segment)
    clear temp;
    temp = zeros(num_procs, num_bounds+1); %create array to store indices for dim i
    proc_count = 1; %keep track of the number of processors
    
    if dim==2 %2D array
        for g1 = 1:len(1) %grid cols
            for g2=1:len(2) %grid rows
                curr_inds = global_ind{g1,g2};
                curr_falls = local_falls{g1,g2};
                if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                    if length([grid(g1,g2) 1 s(dims(i))])==length(temp(proc_count,:))
                        temp(proc_count,:) = [grid(g1,g2) 1 s(dims(i))]; %block range is the entire dimension
                    else %---NEW CODE
                        temp(proc_count, 1:length([grid(g1,g2) 1 s(dims(i))])) = [grid(g1,g2) 1 s(dims(i))];
                    end  %---NEW CODE
                else %dimension is broken up
                    dim_inds = curr_inds{dims(i)};
                    temp2 = [];
                    if curr_falls(dims(i)).complete_cycle
                        for ii=1:curr_falls(dims(i)).n %iterate over the number of line segments in the local falls
                            temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                        end
                    else %incomplete cycle
                        for ii=1:curr_falls(dims(i)).n-1 %iterate over the number of line segments in the local falls
                            %since the block is complete but cycle is not, local
                            %processor has one less line segment
                            temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                        end
                        if ~curr_falls(dims(i)).complete_block
                            %compute last right index
                            last_block_len = rem(curr_falls(dims(i)).local_len, ...
                                (curr_falls(dims(i)).r-curr_falls(dims(i)).l+1));
                            temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1)) ...
                                curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1))+(last_block_len-1)];
                        end
                    end %incomplete cycle
                    %pad with zeros if necessary
                    if length(temp2)<num_bounds
                        num_zeros = num_bounds-(length(temp2));
                        pad = zeros(1,num_zeros);
                        temp2 = [temp2 pad];
                    end
                    temp(proc_count,:) = [grid(g1,g2) temp2]; %!!!!!!
                    clear temp2;
                end %dimension is broken up
                proc_count = proc_count+1;
            end
        end
    elseif dim==3 %3D array
        for g1 = 1:len(1)
            for g2=1:len(2)
                for g3 = 1:len(3)
                    curr_inds = global_ind{g1,g2,g3};
                    curr_falls = local_falls{g1,g2,g3};
                    if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                        temp(proc_count,:) = [grid(g1,g2,g3) 1 s(dims(i))]; %block range is the entire dimension
                    else %block range is not the entire range
                        dim_inds = curr_inds{dims(i)};
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        temp2 = [];
                        if curr_falls(dims(i)).complete_cycle
                            for ii=1:curr_falls(dims(i)).n %iterate over the number of line segments in the local falls
                                temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                    curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                            end
                        else %incomplete cycle
                            for ii=1:curr_falls(dims(i)).n-1 %iterate over the number of line segments in the local falls
                                %since the block is complete but cycle is not, local
                                %processor has one less line segment
                                temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                    curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                            end
                            if ~curr_falls(dims(i)).complete_block
                                %compute last right index
                                last_block_len = rem(curr_falls(dims(i)).local_len, ...
                                    (curr_falls(dims(i)).r-curr_falls(dims(i)).l+1));
                                temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1)) ...
                                    curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1))+(last_block_len-1)];
                            end
                        end %incomplete cycle
                        %pad with zeros if necessary
                        if length(temp2)<num_bounds
                            num_zeros = num_bounds-(length(temp2));
                            pad = zeros(1,num_zeros);
                            temp2 = [temp2 pad];
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        temp(proc_count,:) = [grid(g1,g2,g3) temp2];
                        clear temp2;
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
                        curr_falls = local_falls{g1,g2,g3,g4};
                        if isa(curr_inds{dims(i)}, 'char') & (curr_inds{dims(i)} == ':')
                            temp(proc_count,:) = [grid(g1,g2,g3,g4) 1 s(dims(i))]; %block range is the entire dimension
                        else %block range is not the entire range
                            dim_inds = curr_inds{dims(i)};
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            temp2 = [];
                            if curr_falls(dims(i)).complete_cycle
                                for ii=1:curr_falls(dims(i)).n %iterate over the number of line segments in the local falls
                                    temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                        curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                                end
                            else %incomplete cycle
                                for ii=1:curr_falls(dims(i)).n-1 %iterate over the number of line segments in the local falls
                                    %since the block is complete but cycle is not, local
                                    %processor has one less line segment
                                    temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(ii-1))...
                                        curr_falls(dims(i)).r+(curr_falls(dims(i)).s*(ii-1))];
                                end
                                if ~curr_falls(dims(i)).complete_block
                                    %compute last right index
                                    last_block_len = rem(curr_falls(dims(i)).local_len, ...
                                        (curr_falls(dims(i)).r-curr_falls(dims(i)).l+1));
                                    temp2 = [temp2 curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1)) ...
                                        curr_falls(dims(i)).l+(curr_falls(dims(i)).s*(curr_falls(dims(i)).n-1))+(last_block_len-1)];
                                end
                            end %incomplete cycle
                            %pad with zeros if necessary
                            if length(temp2)<num_bounds
                                num_zeros = num_bounds-(length(temp2));
                                pad = zeros(1,num_zeros);
                                temp2 = [temp2 pad];
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            temp(proc_count,:) = [grid(g1,g2,g3,g4) temp2];
                            clear temp2;
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
