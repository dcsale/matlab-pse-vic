function B = idct(B)
%IDCT Distributed inverse discrete cosine transform.
%   A = IDCT(B) inverts the DCT transform, returning the
%   original vector if B was obtained using B = DCT(A).
%
%   If B is a matrix, the IDCT operation is applied to
%   each column.
%
%   Author: Daniel Jennings (dsj@ll.mit.edu)
%   Date created: June 1, 2005

global pMATLAB;

%make sure input is 2-D
if length(size(B)) ~= 2
    error('dmat/idct.m: Input must be 2-D');
end

%*************************************************
%------> adjust map so it is along columns <------
%*************************************************
%find expected grid spec (want along columns)
new_gridspec = [1 length(B.map.proc_list)];

%if map and new gridspec do not match up then need remap
if (size(B.map.grid, 1) ~= new_gridspec(1)) | ...
   (size(B.map.grid, 2) ~= new_gridspec(2))

   new_proclist = B.map.proc_list;

   %TODO is distspec ok?
   new_distspec(1).dist = 'b';
   new_distspec(2).dist = 'b';

   newMap = map(new_gridspec, new_distspec, new_proclist);
   B = remap(B, newMap);
end
%*************************************************
%----------> done with map adjustment <-----------
%*************************************************

%perform computation
B.local = idct(B.local);


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
