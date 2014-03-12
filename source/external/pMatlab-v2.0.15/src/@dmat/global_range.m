function varargout = global_range(d, varargin)
%GLOBAL_RANGE Returns the ranges of global indices local to the
%   current processor. Returns the same range if D is block distributed,
%   returns subranges for (block) cyclic distirbutions.
%   GLOBAL_RANGE(D, DIM) Returns the global index range of the 
%       distributed array D local to the current processor in the 
%       specified dimension, DIM.
%   GLOBAL_RANGE(D) Returns the global index range of the 
%       distributed array D local to the current processor in 
%       all dimensions of D.
%
% Author:   Nadya Travinin

if length(varargin)==1 %return index range for dimension DIM
    dims = varargin{1};
else %return index range for all dimensions
    dims = 1:length(d.size);
end

my_inds = d.global_ind; %local global indices
my_falls = d.falls;
s = d.size; 
for i = 1:length(dims)
    if isa(my_inds{dims(i)}, 'char') & (my_inds{dims(i)} == ':')
        varargout{i} = [1 s(dims(i))];
    else %dimension is broken up
        dim_inds = my_inds{dims(i)}; %get the indices in the current dimension
        temp = [];
        if my_falls(dims(i)).complete_cycle
            for ii=1:my_falls(dims(i)).n %iterate over the number of line segments in the local falls
                temp = [temp my_falls(dims(i)).l+(my_falls(dims(i)).s*(ii-1))...
                    my_falls(dims(i)).r+(my_falls(dims(i)).s*(ii-1))];
            end
        else %incomplete cycle
            for ii=1:my_falls(dims(i)).n-1 %iterate over the number of line segments in the local falls
                %since the block is complete but cycle is not, local
                %processor has one less line segment
                temp = [temp my_falls(dims(i)).l+(my_falls(dims(i)).s*(ii-1))...
                    my_falls(dims(i)).r+(my_falls(dims(i)).s*(ii-1))];
            end
            if ~my_falls(dims(i)).complete_block
                %compute last right index
                last_block_len = rem(my_falls(dims(i)).local_len, ...
                    (my_falls(dims(i)).r-my_falls(dims(i)).l+1));
                temp = [temp my_falls(dims(i)).l+(my_falls(dims(i)).s*(my_falls(dims(i)).n-1)) ...
                    my_falls(dims(i)).l+(my_falls(dims(i)).s*(my_falls(dims(i)).n-1))+(last_block_len-1)];
            end
        end %incomplete cycle
            
        varargout{i} = temp;
        clear temp;
    end %dimension is broken up
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
