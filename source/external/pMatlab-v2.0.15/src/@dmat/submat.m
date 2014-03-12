function b = submat(a, s)
% SUBMAT  Helper function used by subsref.
%   SUBMAT(A, S) Helper function used by subscripted reference. A is a
%   distributed array, S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts.
%
% Author:   Nadya Travinin
%

global vars
global pMATLAB;

if a.dim==2  %2-D
   if isa(s.subs{1}, 'char') & isa(s.subs{2}, 'char')
       if (s.subs{1}==':') & (s.subs{2}==':') %A(:,:)
           b = a;
       end
   elseif (length(s.subs{1})==1) & (length(s.subs{2})==1) ...
           & (~isa(s.subs{1}, 'char')) & (~isa(s.subs{2}, 'char')) %A(i,j)
       %----CASE 1 on SUBSREF development plan----
       ind = get_ind_range(a,s);
       local_ind = get_local_ind(a.global_ind, ind);
       
       %need to figure out on all procs where the index is local
       r = get_local_proc(a.pitfalls, a.map.grid, [ind{1} ind{2}]);
       %define new map
       new_map = map([1 1], {}, [r]);
       
       %create a new dmat
       b.map = new_map;
       b.dim = a.dim;
       b.size = [length(ind{1}) length(ind{2})];
       %create a PITFALLS for each dimension
       for i=1:new_map.dim
           b.pitfalls(i) = gen_pitfalls(size(new_map.grid,i), new_map.dist_spec(i), 1);
       end
       b.falls = get_local_falls(b.pitfalls, new_map.grid, pMATLAB.my_rank);
       b.local = a.local(local_ind{1}, local_ind{2});
       grid_dims = size(new_map.grid);
       b.global_ind = get_global_ind(b.falls, grid_dims);
       
       b = class(b, 'dmat');
       %----CASE 1 on SUBSREF development plan----
   else
       ind = get_ind_range(a,s);
       local_ind = get_local_ind(a.global_ind, ind);
       b.map = a.map;
       b.dim = a.dim;
       b.size = [length(ind{1}) length(ind{2})];
       b.pitfalls = a.pitfalls;
       b.falls = a.falls;
       b.local = a.local(local_ind{1}, local_ind{2});
       b.global_ind = a.global_ind;
       
       b = class(b, 'dmat');    
   end
elseif a.dim==3  %3-D
   if (s.subs{1}==':')&(s.subs{2}==':')&(s.subs{3}==':') %A(:,:,:)
       b = a;
   else
       ind = get_ind_range(a,s);
       local_ind = get_local_ind(a.global_ind, ind);
       b.map = a.map;
       b.dim = a.dim;
       b.size = [length(ind{1}) length(ind{2}) length(ind{3})];
       b.pitfalls = a.pitfalls;
       b.falls = a.falls;
       b.local = a.local(local_ind{1}, local_ind{2}, local_ind{3});
       b.global_ind = a.global_ind;
       
       b = class(b, 'dmat');
   end
elseif a.dim==4 %4-D
    if (s.subs{1}==':')&(s.subs{2}==':')&(s.subs{3}==':')&(s.subs{4}==':') %A(:,:,:,:)
       b = a;
   else
       ind = get_ind_range(a,s);
       local_ind = get_local_ind(a.global_ind, ind);
       b.map = a.map;
       b.dim = a.dim;
       b.size = [length(ind{1}) length(ind{2}) length(ind{3}) length(ind{4})];
       b.pitfalls = a.pitfalls;
       b.falls = a.falls;
       b.local = a.local(local_ind{1}, local_ind{2}, local_ind{3}, local_ind{4});
       b.global_ind = a.global_ind;
       
       b = class(b, 'dmat');
   end
else
   error('SUBMAT2: Only up to 4 dimensions are supported.');
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
