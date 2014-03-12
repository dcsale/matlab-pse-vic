function p = gen_pitfalls(np, dist_spec, dim_len, varargin)
%GEN_PITFALLS Given the number of processors, distribution spec, and the
%   length of the dimension, generates all the PITFALLS information.
%   GEN_PITFALLS(NP, DIST_SPEC, DIM_LEN)
%       NP - number of processors along which the dimension is distributed
%       DIST_SPEC - distribution specification with two possible fields.
%           The DIST field is mandatory and could have the following values:
%           'b' - block
%           'c' - cyclic
%           'bc' - block-cyclic
%           If DIST_SPEC.DIST == 'bc', the block size 'B_SIZE' must
%           also be defined.
%       DIM_LEN - length of the dimension for which the PITFALLS is being
%       created
%
%   Returns p:
%       P is a single PITFALLS used to describe the distribution. P will
%       have an extra field (P.REM_CYCLE) which will store the length of
%       the incomplete cycle in the current distribution. This way 3
%       PITFALLS will not be necessary to represent non-evenly divisible
%       dimensions.
%       P has the following fields:
%           p.L - starting index of the first block on the first processor
%           p.R - ending index of the first block on the first processor
%           p.S - stride between successive L's
%           p.N - number of equally spaced, equally sized blocks of elements
%               per processor
%           p.D - spacing between L's of successive processor FALLS
%           p.P - number of processors
%           p.REM_CYCLE - the length of the incomplete cycle. IMPORTANT:
%               p.N includes the incomplete cycle if one exists, i.e. the
%               PITFALLS calculation rounds up the length of the dimension and
%               computes the PITFALLS information as if the dimension was
%               evenly divisible by the cycle length.
%
% Author:   Nadya Travinin
%
% References: Shankar Ramaswamy and Prithviraj Banerjee. Automatic Generation of Efficient Array 
%             Redistribution Routines for Distributed Memory Multicomputers. IEEE 1995.

if length(varargin) == 0  %no overlap
    %store block size
    if strcmp(dist_spec.dist, 'b')
        b_size = ceil(dim_len/np);
    elseif strcmp(dist_spec.dist, 'c')
        b_size = 1;
    elseif strcmp(dist_spec.dist, 'bc')
        b_size = dist_spec.b_size;
    end
    
    %cycle length
    cycle_len = b_size*np;
    %number of cycles - both complete and incomplete
    num_cycles = ceil(dim_len/cycle_len);
    %length of incomplete cycle
    rem_cycle = mod(dim_len, cycle_len);
    
    %create the PITFALLS data structure
    p.l = 1;
    p.r = b_size;
    
    %!!!THIS PRODUCES STRIDE>1 EVEN IF NUM_CYCLES==1...MAKE SURE THIS WORKS
    %!!!This might influence the intersection algorithm
    p.s = b_size*np; 
    
    p.n = num_cycles;
    p.d = b_size;
    p.p = np;
    p.rem_cycle = rem_cycle;    
else  %overlap
    overlap = varargin{1};
    %store block size
    if strcmp(dist_spec.dist, 'b')
        b_size = ceil(dim_len/np);
    else
        error('GEN_PITFALLS: Overlap is only supported for block distributions.');
    end
    
    %cycle length
    cycle_len = b_size*np+overlap; 
    %number of cycles - both complete and incomplete
    num_cycles = ceil(dim_len/cycle_len);
    %length of incomplete cycle
    rem_cycle = mod(dim_len, cycle_len);
    
    %create the PITFALLS data structure
    p.l = 1;
    p.r = b_size+overlap;
    
    %!!!THIS PRODUCES STRIDE>1 EVEN IF NUM_CYCLES==1...MAKE SURE THIS WORKS
    %!!!This might influence the intersection algorithm
    p.s = b_size*np; 
    
    p.n = num_cycles;
    p.d = b_size;
    p.p = np;
    p.rem_cycle = rem_cycle;    
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