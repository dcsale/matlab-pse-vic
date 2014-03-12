function f = pitfalls_intersection(p1,p2)
%PITFALLS_INTERSECTION Computes intersection of two PITFALLS.
%   P1 and P2 can consist of multiple PITFALLS since a maximum of three
%   PITFALLS are needed to describe a distribution. F is an 2D array of
%   structures, where F(i,j) contains information about what communication
%   needs to occur between i-th processor of P1 and j-th processors of P2.
%   That information is stored in a structure C with the following fields:
%       C.FALLS - FALLS data structure (or array of FALLS data structures
%           with fields (L, R, S, N)
%       C.PROCS - array of indices of the two processors that need to communicate data
%           indices defined by C.FALLS.
%
% Author:   Nadya Travinin
%
% References: Shankar Ramaswamy and Prithviraj Banerjee. Automatic Generation of Efficient Array 
%             Redistribution Routines for Distributed Memory Multicomputers. IEEE 1995.


%compute the number of PITFALLS used to describe each of the two
%distributions
len1 = length(p1);
len2 = length(p2);

if (len1==1) & (len2==1)
    for i1 = 0:p1.p-1
        f1.l = p1.l+i1*p1.d;
        f1.r = p1.r+i1*p1.d;
        f1.s = p1.s;
        f1.n = p1.n;
        for i2 = 0:p2.p-1
            f2.l = p2.l+i2*p2.d;
            f2.r = p2.r+i2*p2.d;
            f2.s = p2.s;
            f2.n = p2.n;
            fi = falls_intersection(f1,f2);
            if ~isempty(fi)
                c.falls = fi;
                %store (grid) INDICES of communicating processors in the
                %current dimension
                c.procs = [i1+1 i2+1];
                f(i1+1, i2+1) =  c;
            end
        end
    end
else
    %!!!FOR NOW - only implement the algorithm for FALLS objects that consist
    %of a single FALLS (i.e. evenly divisible dimensions)
    error('PITFALLS_INTERSECTION: PITFALLS objects only of length 1 currently supported');
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