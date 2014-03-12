function local_size = localdims(falls, dim)
%LOCALDIMS Given the local FALLS and object dimension, calculates
%   a vector of local dimensions of the distributed object.
%   LOCALDIMS(FALLS, DIM)
%       FALLS - an array of FALLS structures for each dimension of
%       the distributed object, i.e. FALLS(i) is a FALLS
%       representation of the i-th dimension
%       DIM - dimension of the distributed object (2, 3, or 4)
%
%   Returns LOCAL_SIZE:
%       length(LOCAL_SIZE)is equal to the number of dimensions of the
%       distributed object.
%
% Author:   Nadya Travinin

if dim<=4
    local_falls = falls;
    %compute local dimensions from the local FALLS
    %!!!NOTE: The following calculation works only if data on processor
    %RANK can be represented by a single FALLS.
    local_size = [];
    if isa(local_falls(1), 'struct')
        for i = 1:dim
            %local size is just the local_len field of the local falls
            local_size = [local_size local_falls(i).local_len];
        end
    else %no local data
        local_size = zeros(1, dim);
    end
else
    error('LOCALDIMS: Only objects up to 4-D are supported');
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