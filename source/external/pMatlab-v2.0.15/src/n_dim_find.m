function ind = n_dim_find(array, val)
%N_DIM_FIND N-dimensional FIND function that returns array of indices for
%   val in each dimensions. IND is empty matrix [] if VAL was not
%   found in ARRAY
%
% Author:   Nadya Travinin

%!!!Only supports 4-D arrays at the moment

%dimension of the array
dim = length(size(array));

if dim>4
    error('N_DIM_FIND: Only four dimensional objects supported');
end

[i, j] = find(array==val);
if isempty(i) && isempty(j)
    ind = [];
elseif dim==2
    ind(1) = i;
    ind(2) = j;
elseif dim==3
    ind(1) = i;
    num_cols = size(array,2);
    temp = mod(j, num_cols);
    if temp==0
        ind(2) = num_cols;
    else
        ind(2) = temp;
    end
    ind(3) = ceil(j/num_cols);
elseif dim==4
    ind(1) = i;
    num_cols = size(array,2);
    temp = mod(j, num_cols);
    if temp==0
         ind(2) = num_cols;
    else
        ind(2) = temp;
    end
    num_slices = ceil(j/num_cols);
    dim3 = size(array,3);
    temp = mod(num_slices, dim3);
    if temp==0
        ind(3) = dim3;
    else
        ind(3) = temp;
    end
    ind(4) = ceil(num_slices/dim3);
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
