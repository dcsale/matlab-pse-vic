function local_falls = get_local_falls(pitfalls, grid, rank)
%GET_LOCAL_FALLS Given the PITFALLS object, the grid and local processor
%   rank, returns an array of local FALLS objects (one for each dimension of the
%   distributed object).
%   GET_LOCAL_FALLS(PITFALLS, GRID, RANK)
%       PITFALLS - an array of PITFALLS structures for each dimension of
%       the distributed object, i.e. PITFALLS(i) is a PITFALLS
%       representation of the i-th dimension
%       GRID - n-dimensional array that specifies the layout of processors
%       on which the object is distributed
%       RANK - processor rank
%
%   Returns LOCAL_FALLS:
%       Array of FALLS objects, one for each dimension of the distributed
%       object. FALLS(i) is a 4-tuple (l,r,s,n). In addition FALLS(i)
%       stores information regarding local length of data to take care of
%       cases with non-evenly divisible dimensions without implementing
%       extra FALLS. 
%       The following are the fields of each FALLS(i)
%           L - starting index of the first block
%           R - ending index of the first block
%           S - stride between successive L's
%           N - number of equally spaced, equally sized blocks of elements
%           LOCAL_LEN - local length of the data on the local processor
%               NOTE: If dimensions are evenly divisible, LOCAL_LEN is
%               always (R-L+1)*N
%           COMPLETE_CYCLE - flag that indicates whether the processor has
%           incomplete cycles
%           COMPLETE_BLOCK - flag that indicates whether the processor has
%           incomplete blocks
%
% Author:   Nadya Travinin
%
% References: Shankar Ramaswamy and Prithviraj Banerjee. Automatic Generation of Efficient Array 
%             Redistribution Routines for Distributed Memory Multicomputers. IEEE 1995.

dim = length(pitfalls);

if dim <= 4
    ind = n_dim_find(grid, rank);
    if (~isempty(ind)) && (length(ind)<dim)
        for d = (length(ind)+1):dim
            ind(d) = 1;
        end
    end
    if ~isempty(ind)
        for i = 1:dim
            f.l = pitfalls(i).l+(ind(i)-1)*pitfalls(i).d;
            f.r = pitfalls(i).r+(ind(i)-1)*pitfalls(i).d;
            f.s = pitfalls(i).s;
            f.n = pitfalls(i).n;
            
            %compute local length
            block_size = f.r-f.l+1; %block size 
            
            if block_size == pitfalls(i).d %no overlap
                if pitfalls(i).rem_cycle==0
                    f.local_len = f.n*block_size;
                    f.complete_cycle = logical(1); %flag that signifies that the local data has no incomplete cycles
                    f.complete_block = logical(1); %flag that signifies that the local data has no incomplete blocks
                    %(i.e. local indices can be computed directly from the falls info without 
                    % dealing with local length)
                else %rem_cycle ~= 0
                    %number of blocks in the incomplete cycle (both complete
                    %and incomplete)
                    num_blocks = ceil(pitfalls(i).rem_cycle/block_size);
                    %the length of incomplete block (if one exists)
                    rem_block = mod(pitfalls(i).rem_cycle, block_size); 
                    if ind(i)<num_blocks  
                        %local processor does not have an incomplete cycle
                        %and has no incomplete blocks
                        f.local_len = f.n*block_size;
                        f.complete_cycle = logical(1);
                        f.complete_block = logical(1);
                    elseif ind(i)>num_blocks
                        %local processor has one less block (incomplete cycle)
                        %but no imcomplete blocks
                        f.local_len = (f.n-1)*block_size;  
                        f.complete_cycle = logical(0);
                        f.complete_block = logical(1);
                    elseif ind(i)==num_blocks
                        if rem_block==0
                            %local processor has no incomplete cycles or blocks
                            f.local_len = f.n*block_size;
                            f.complete_cycle = logical(1);
                            f.complete_block = logical(1);
                        else
                            %local processor has an incomplete cycle with an
                            %incomplete block
                            f.local_len = (f.n-1)*block_size+rem_block;
                            f.complete_cycle = logical(0);
                            f.complete_block = logical(0);
                        end
                    end
                end
            else %overlap, rem_cycle should always be greater than 0 
                if f.r <= pitfalls(i).rem_cycle
                    f.local_len = block_size;
                    f.complete_cycle = logical(1);
                    f.complete_block = logical(1);
                else
                    f.local_len = pitfalls(i).rem_cycle - f.l+1;
                    f.complete_cycle = logical(0);
                    f.complete_block = logical(0);
                end
            end %overlap
            
            %store the FALLS structure 
            local_falls(i) = f;
            
            clear f;
        end
    else
        for i = 1:dim
            local_falls(i) = -1;
        end
    end
else
    error('get_local_falls: Only objects up to 4-D are supported');
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
