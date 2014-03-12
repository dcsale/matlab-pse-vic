function local_ind = get_local_ind(global_ind, ind)
%function local_ind = get_local_ind(grid, global_ind, ind)
%GET_LOCAL_IND Returns a cell array of local indices given an array of
%   global indices on the current processor and an array of global indices
%   that are being referenced. 
%   GET_LOCAL_IND(GLOBAL_IND, IND)
%       GLOBAL_IND - array of length equal to the number of dimensions.
%       Each entry in the array specifies the global indices in the i-th
%       dimension stored on the current processor.
%       IND - array of length equal to the number of dimensions, where each
%       entry specifies the global indices being referenced in that
%       dimension.
%           
%   Returns LOCAL_IND:
%       length(LOCAL_IND) is equal to the number of dimensions of the distributed
%       object. LOCAL_IND(i) is of the form [ind1 ind2 ind3 ...] where ind_i is a
%       local index of the data stored locally.
%
% Author:   Nadya Travinin

dim = length(ind);

if dim <= 4 %number of dimensions of the distributed object is 4 or less
    for i = 1:dim
        local_ind{i} = [];
        if strcmp(global_ind{i}, ':') %dimension i is not distributed
            loc_inds = ind{i};
        else  %dimension i is distributed
            if strcmp(ind{i}, ':')
                loc_inds = ':';
            else
                %PRESERVES THE ORDERING OF INDICES
                 loc_inds = [];
%                 for j=1:length(ind{i})
%                     temp_ind = find(global_ind{i}==ind{i}(j));
%                     if ~isempty(temp_ind)
%                         loc_inds = [loc_inds temp_ind];
%                     end
%                 end

                % vvvvv new code vvvvv
                [vals i_global i_ind] = intersect(global_ind{i}, ind{i});
                [ind_sorted i_ind_sorted] = sort(i_ind);
                loc_inds = i_global(i_ind_sorted);
                if (length(loc_inds) == 0)
                   loc_inds = [];
                end
                % ^^^^^ new code ^^^^^

              end
        end   %dimension i is distributed
        local_ind{i} = loc_inds;
    end
else  %number of dimensions is greater than 4
    error('localdims: Only objects up to 4-D are supported');
end   %number of dimensions is greater than 4

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