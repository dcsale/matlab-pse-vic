function d = dmat(varargin)
% DMAT Distributed matrix constructor.
%   Creates the necessary data structures for the distributed matrix.
%
%   NOTE: DMAT does not allocate memory for the distributed matrix
%   and should never be called directly.  A distributed matrix should
%   always be created with an appropriate constructor (ZEROS, ONES, RAND
%   or SPARSE) which will be responsible for allocating memory.
%
%   DMAT(N, P) If N is a scalar, creates the necessary data
%       structures for an NxN distributed matrix; if N is a vector, 
%       creates the necessary structures for an N distributed matrix
%   DMAT(M, N, P) Creates the necessary data structures for an MxN
%       distributed matrix.
%   DMAT(M, N, Q, P) Creates the necessary data structures for an MxNxQ
%       distributed array.
%   DMAT(M, N, Q, R, P) Creates the necessary data structures for an
%       MxNxQxR distributed array.
%   
%   Returns D:
%   D.MAP - processor map onto which the distributed array is mapped.  
%   D.DIM - dimension of the distributed array. The dimension of the map has 
%       to be consistent with the dimension of the distributed array.
%   D.SIZE - size of the distributed array, i.e. number of elements in each
%       dimension.
%   D.PITFALLS - pitfalls structure that describes the distribution of the
%       distributed array among the processors in the map.
%   D.FALLS - falls structure that describes the distribution of the
%       distributed object on the local processor.
%   D.LOCAL - numerical part of the distributed array stored on the local
%       processor.  Memory for D.LOCAL will be allocated by constructor
%       functions, such as ZEROS, ONES, RAND and SPARSE.
%   D.GLOBAL_IND - global indices of the distributed array local to the
%       current processor.
%
%   See also: ONES, ZEROS, RAND, SPARSE.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

superiorto('double');

%global vars
global pMATLAB;

dims = [varargin{1:end-1}];
p = varargin{end};

if length(dims) == 0
  error('@dmat/dmat: Must specify at least one dimension');
elseif length(dims) == 1 %DMAT(M, P)
  dims = [dims dims];
elseif length(dims) >= 5 
  error('@dmat/dmat: Incorrect number of inputs');  
end

d.map = p;
d.dim = length(dims);
d.size = dims;
if p.dim ~= length(dims)
  error('@dmat/dmat: Map and distributed object dimensions must match');
end

%create a PITFALLS for each dimension
for i=1:p.dim
  if isempty(p.overlap)
    d.pitfalls(i) = gen_pitfalls(size(p.grid,i), p.dist_spec(i), dims(i));
  elseif p.overlap(i)==0
    d.pitfalls(i) = gen_pitfalls(size(p.grid,i), p.dist_spec(i), dims(i));
  else
    d.pitfalls(i) = gen_pitfalls(size(p.grid,i), p.dist_spec(i), dims(i), p.overlap(i));
  end
end
  
%get the local falls
d.falls = get_local_falls(d.pitfalls, p.grid, pMATLAB.my_rank);
  
%figure out local dimensions
local_size = localdims(d.falls, d.dim);
  
% Allocating memory is the responsibility of map functions
% (e.g. ones, zeros, rand and sparse)
%    d.local = zeros(local_size);
d.local = [];
  
%get the local indices for the current processor
grid_dims = size(p.grid);
if length(grid_dims)<p.dim
  for i = (length(grid_dims)+1):p.dim
    grid_dims(i) = 1;
  end
end
d.global_ind = get_global_ind(d.falls, grid_dims);
d = class(d, 'dmat'); 

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
