function fi = falls_intersection(f1, f2)
%FALLS_INTERSECTION Given two FALLS, find the intersection FALLS
%   F1 and F2 each contain at most a single FALLS, since unevenly
%   divisible dimensions are represented using a single FALLS and local
%   length.
%
%   FALLS data structure F:
%       F.L - beginning index of the first block of elements
%       F.R - ending index of the first block of elements
%       F.S - strides between successive L's
%       F.N - number of equally spaced, equally sized blocks
%
% Author:   Nadya Travinin
%
% References: Shankar Ramaswamy and Prithviraj Banerjee. Automatic Generation of Efficient Array 
%             Redistribution Routines for Distributed Memory Multicomputers. IEEE 1995.

% check to see if either f1 or f2 are -1's
if ~isstruct(f1) || ~isstruct(f2)
  fi = [];
  return
end

%compute intersection period (fp), m1 and m2
fp = lcm(f1.s, f2.s);
m1 = fp/f1.s;
m2 = fp/f2.s;

fi = [];
I1 = max(0, ceil((f2.l-f1.r)/f1.s));
for i1=I1 : min(I1+m1-1, f1.n-1)
    l1.l = f1.l+i1*f1.s;
    l1.r = f1.r+i1*f1.s;
    for i2 = max(0, ceil((i1*f1.s+f1.l-f2.r)/f2.s)) : min([floor((i1*f1.s+f1.r-f2.l)/f2.s) m2-1 f2.n-1])
        l2.l = f2.l+i2*f2.s;
        l2.r = f2.r+i2*f2.s;
        li = ls_intersection(l1,l2);
        %don't have to check that intersection is non-empty since we
        %are guaranteed intersection by iterating over the above
        %specified loop bounds
        temp.l = li.l;
        temp.r = li.r;
        temp.s = fp;
        temp.n = floor(min((f1.n-i1-1)/m1, (f2.n-i2-1)/m2) + 1);
        %last global index of the intersection
        fi_block_size = temp.r-temp.l+1; %block size
        fi_last_ind = temp.r+temp.s*(temp.n-1); 
        %last global index of f1
        if (f1.complete_block & f1.complete_cycle) 
            f1_last_ind  = f1.r+f1.s*(f1.n-1);
        elseif f1.complete_block
            f1_last_ind = f1.r+f1.s*(f1.n-2);
        else
            f1_block_size = f1.r-f1.l+1;
            f1_rem_block = mod(f1.local_len, f1_block_size);
            f1_last_ind = f1.r + f1.s*(f1.n-1) - (f1_block_size-f1_rem_block);
        end
        %last global index of f2
        if (f2.complete_block & f2.complete_cycle) 
            f2_last_ind  = f2.r+f2.s*(f2.n-1);
        elseif f2.complete_block
            f2_last_ind = f2.r+f2.s*(f2.n-2);
        else
            f2_block_size = f2.r-f2.l+1;
            f2_rem_block = mod(f2.local_len, f2_block_size);
            f2_last_ind = f2.r+f2.s*(f2.n-1) - (f2_block_size-f2_rem_block);
        end
        
        if (fi_last_ind <= f1_last_ind) & (fi_last_ind <= f2_last_ind)
            temp.local_len = temp.n*fi_block_size;
            temp.complete_cycle = logical(1);
            temp.complete_block = logical(1);
        else
            last_ind = min(f1_last_ind, f2_last_ind);
            diff_ind = fi_last_ind - last_ind;
            temp.local_len = temp.n*fi_block_size - diff_ind;
            if diff_ind >= fi_block_size %!!!
                temp.complete_cycle = logical(0);
                temp.complete_block = logical(1);
                
            elseif diff_ind < fi_block_size
                temp.complete_cycle = logical(0);
                temp.complete_block = logical(0);
            end
            
        end
        
        fi = [fi temp];
    end    
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