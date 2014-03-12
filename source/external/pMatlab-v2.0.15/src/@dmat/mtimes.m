function c = mtimes(a,b)
% *  Matrix multiply.
%   C = MTIMES(A, B) or C = A*B multiples two matrices together. 
%   Mtimes is an inner product matrix multiply and assumes that A is row mapped and B is column mapped.  
%   If B is row mapped, it is remapped in order to perform the matrix multiplication.   
%
%   For A row mapped, the function accepts B matrices that are double (not mapped), 
%    column mapped or row mapped.  
%   For a matrix A that is not mapped (doubles), B must be column mapped. 
%
%   WARNING: Overlaps have not been tested.  In fact, distributed
%   matrix multiply was implemented without a thought to overlap.
%   This is a major TODO.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

global pMATLAB;

% If one of the arguments is a scalar, do element-wise multiplication
if (prod(size(a)) == 1 | prod(size(b)) == 1)
   c = a .* b;
   return;
end

% Check that both arguments are 2D
if (length(size(a)) ~= 2) | (length(size(b)) ~= 2)
   error('DMAT/MTIMES: Input arguments must be 2-D');
end

if isa(a, 'double') && isa(b, 'dmat')
  % double * dmat
  mapB = b.map;
  
  % create sub-result matrices, each with a map equivalent to 
  % one of map B's rows
  res_gridspec = [1 size(mapB.grid, 2)];
  res_distspec(1).dist = 'b';
  res_distspec(2) = mapB.dist_spec(2);
  my_idx = 1;
  for i=1:size(mapB.grid, 1)
    res_map = map(res_gridspec, res_distspec, mapB.grid(i, :));
    res{i} = zeros(size(a, 1), size(b, 2), res_map);

    % store which result matrix this processor will use
    if inmap(res_map, pMATLAB.my_rank)
      my_idx = i;
    end
  end  
  
  % map A's ith column to the leftmost node in B's ith row
  new_gridspec = [1 size(mapB.grid, 1)];
  new_distspec(1).dist = 'b';
  new_distspec(2) = mapB.dist_spec(1);
  new_proclist = mapB.grid(:, 1)';
  if ~isempty(mapB.overlap)
    warning('@dmat/mtimes: TODO dmats with overlap may not work');
    new_mapA = map(new_gridspec, new_distspec, new_proclist, mapB.overlap);
  else
    new_mapA = map(new_gridspec, new_distspec, new_proclist);
  end
  a = remap(a, new_mapA);

  % find row/column in grid
  [myRow, myCol] = find(mapB.grid == pMATLAB.my_rank);
    
  % send/receive data and do the multiplication
  res{my_idx}.local = multicast(mapB.grid(myRow, 1), ...
                                mapB.grid(myRow, 2:end), ...
                                a.local) * b.local;

  % add back together sub-results to form resulting matrix (c)
  c = summation(zeros(size(a, 1), size(b, 2), mapB), res);
  
else
  % dmat * dmat OR dmat * double
  mapA = a.map;

  % create sub-result matrices, each with a map equivalent to 
  % one of map A's columns
  res_gridspec = [size(mapA.grid, 1) 1];
  res_distspec(1) = mapA.dist_spec(1);
  res_distspec(2).dist = 'b';
  my_idx = 1;
  for i=1:size(mapA.grid, 2)
    res_map = map(res_gridspec, res_distspec, mapA.grid(:, i));
    res{i} = zeros(size(a, 1), size(b, 2), res_map);
    % store which result matrix this processor will use
    if inmap(res_map, pMATLAB.my_rank)
      my_idx = i;
    end
  end
  
  % map B's ith row to the top node in A's ith column
  new_gridspec = [size(mapA.grid, 2) 1];
  new_distspec(1) = mapA.dist_spec(2);
  new_distspec(2).dist = 'b';
  new_proclist = mapA.grid(1, :)';
  if ~isempty(mapA.overlap)
    warning('@dmat/mtimes: TODO dmats with overlap may not work');
    new_mapB = map(new_gridspec, new_distspec, new_proclist, mapA.overlap);
  else
    new_mapB = map(new_gridspec, new_distspec, new_proclist);
  end
  b = remap(b, new_mapB);

  % find row/column in grid
  [myRow, myCol] = find(mapA.grid == pMATLAB.my_rank);
    
  % send/receive data and do the multiplication
  res{my_idx}.local = a.local * multicast(mapA.grid(1, myCol), ...
                                          mapA.grid(2:end, myCol), ...
                                          b.local);

  % add back together sub-results to form resulting matrix (c)
  c = summation(zeros(size(a, 1), size(b, 2), mapA), res);
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
