function b = transpose(x)
%.' Transpose.
%    X.' is the non-conjugate transpose of a distributed matrix. The input
%    ditributed array must be a matrix. The distribution of the distirbuted
%    array is limited to block in 1, 2nd, or both dimensitons.
%
% Author: Nadya Travinin

if (x.dim==2) %MATRIX TRANSPOSE
    xmap = x.map;
    %check that the underslying ditribution is BLOCK, otherwise throw an
    %error
    xdist_spec = xmap.dist_spec;
    for ii = 1:2
        if xdist_spec(ii).dist ~= 'b'
            error('dmat/transpose: Transpose operation currently supported for BLOCK distirbution only.')
        end
    end
    %transpose the map
    tgrid = xmap.grid';
    tproc_list = tgrid(:);
    tmap = map(size(xmap.grid'), {}, tproc_list);
    %create a temp matrix with transpose map
    xsize = x.size;
    %corner turn the original matrix
    xtemp = remap(x, tmap);
    
    %create a second temporary matrix to store the result
    xtemp_out = zeros(xsize(2), xsize(1), xmap);
    
    %store the result (transpose the local part)
    xtemp_out.local(:,:) = (xtemp.local).';
    %return the transpose
    b = xtemp_out;
else
    error('DMAT/TRANSPOSE: Can only be called on matrices. See PERMUTE for higher dimensions.');
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