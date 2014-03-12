function ind = get_global_ind(falls, varargin)
%GET_GLOBAL_IND Returns a cell array of global indices stored locally given
%   a FALLS object (array of FALLS, one for each dimension)
%   GET_GLOBAL_IND(FALLS)
%       FALLS - array of FALLS objects
%       
%   Returns IND:
%       length(IND) is equal to the number of dimensions of the distributed
%       object. IND(i) is of the form [ind1 ind2 ind3 ...] where ind_i is a
%       global index of the distributed object that is local to the current
%       processor.
%
% Author:   Nadya Travinin

dim = length(falls);
if nargin==2
    grid_dims = varargin{1};
else
    grid_dims=[];
end

if dim<=4
    for i = 1:dim
        temp = [];
        if isa(falls(i), 'struct')
            if (~isempty(grid_dims)) && (grid_dims(i) == 1) %dimension is not distributed
                temp = ':';
            else %dimension is distributed
                %get the indices for the first n-1 cycles
                %NOTE: This takes care of the case when complete_cycle==0 and
                %      complete_block==1
                for j=1:falls(i).n-1 
                    temp = [temp falls(i).l+(j-1)*falls(i).s : falls(i).r+(j-1)*falls(i).s];
                end
                %get indices for the n-th cycle, if it exists
                if (falls(i).complete_cycle && falls(i).complete_block) %complete n-th cycle, complete block
                    temp = [temp falls(i).l+(falls(i).n-1)*falls(i).s : falls(i).r+(falls(i).n-1)*falls(i).s];   
                elseif (~falls(i).complete_cycle && ~falls(i).complete_block)
                    %incomplete n-th cycle, incomplete block 
                    block_size = falls(i).r-falls(i).l+1;
                    rem_block = mod(falls(i).local_len, block_size);
                    temp = [temp falls(i).l+(falls(i).n-1)*falls(i).s : falls(i).l+(falls(i).n-1)*falls(i).s + rem_block-1];
                end
            end %dimension is distributed
        end
        ind{i}=temp;
    end
else
    error('GET_GLOBAL_IND: Only objects up to 4-D are supported');
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
