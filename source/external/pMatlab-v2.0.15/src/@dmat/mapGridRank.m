function [gridIndex] = mapGridRank(d)
% Identify grid map index for a given rank and returns gridIndex array
%   Pid ranks:    0 -> Np-1
%   Pid posiiton: 1 -> Np
%
% Author: Dr. Chansup Byun
% Date:   May 3, 2010

% Determine the size of each dimension
dim = size(d.map.grid);
if length(dim)<d.dim
   for i = (length(dim)+1):d.dim
       dim(i) = 1;
   end
end

gridIndex = zeros(Np,d.dim);
if d.dim == 2
   for i = 1:dim(1)
       for j = 1:dim(2)
           gridIndex(d.map.grid(i,j)+1,1) = i;
           gridIndex(d.map.grid(i,j)+1,2) = j;
       end
   end
elseif d.dim==3
   for i = 1:dim(1)
       for j = 1:dim(2)
           for k = 1:dim(3)
               gridIndex(d.map.grid(i,j,k)+1,1) = i;
               gridIndex(d.map.grid(i,j,k)+1,2) = j;
               gridIndex(d.map.grid(i,j,k)+1,3) = k;
           end
       end
   end
elseif d.dim==4
   for i = 1:dim(1)
       for j = 1:dim(2)
           for k = 1:dim(3)
               for m = 1:dim(4)
                   gridIndex(d.map.grid(i,j,k,m)+1,1) = i;
                   gridIndex(d.map.grid(i,j,k,m)+1,2) = j;
                   gridIndex(d.map.grid(i,j,k,m)+1,3) = k;
                   gridIndex(d.map.grid(i,j,k,m)+1,4) = m;
               end
           end
       end
   end
else
   error('DMAT/mapGridRank: Only up to 4-D objects currently supported');
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
