function A = transpose_grid(B);
%TRANSPOSE_GRID Redistributes a dmat by transposing its grid.
%   A = TRANSPOSE_GRID(B) creates a dmat A that has the same contents as the
%   dmat B, except that the grid for A's map is the transpose of the grid for
%   B's map.  The contents of B are automatically redistributed to A.
%   TRANSPOSE_GRID is only supported for 2D dmats.
%
%   TRANSPOSE_GRID is optimized to redistribute row-distributed dmats into
%   column-distributions and vice versa.  For example, suppose B's map had a
%   grid specification of [1 4].  Then A's map will have a grid specification
%   of [4 1].  Note that B must be block distributed in both dimensions with
%   no overlap.
%
%   For all other distributions, e.g. dmats with row and column
%   distributions, overlap, etc., SUBSASGN will be used to redistribute B.
%

global pMATLAB;

comm    = pMATLAB.comm;
my_rank = pMATLAB.my_rank;

% Check that B is 2D
if (ndims(B) ~= 2)
   error('TRANSPOSE_GRID for ND array is not defined.');
end

% If B is not a dmat, simply do a copy
if (~strcmp(class(B), 'dmat'))
   A = B;

% If B is a dmat, perform the redistribution
else

   % Get B's size, map, grid, dist, overlap and cpus list.
   Bsize      = size(B);
   Bmap       = B.map;
   Bgrid      = Bmap.grid;
   Bgrid_spec = size(Bgrid);
   Bdist_spec = Bmap.dist_spec;
   Bproc_list = Bmap.proc_list;
   Boverlap   = Bmap.overlap;

   % Create map for A that is the same as B's but with a transposed grid.
   Agrid_spec = [Bgrid_spec(2) Bgrid_spec(1)];
   if (isempty(Boverlap))
      mapA = map(Agrid_spec, Bdist_spec, Bproc_list);
   else
      mapA = map(Agrid_spec, Bdist_spec, Bproc_list, Boverlap);
   end

   % Construct A
   switch (length(Bsize))
      case {2}
         A = zeros(Bsize(1), Bsize(2), mapA);

      case {3}
         A = zeros(Bsize(1), Bsize(2), Bsize(3), mapA);

      case {4}
         A = zeros(Bsize(1), Bsize(2), Bsize(3), Bsize(4), mapA);

   end

   % Pick off everthing but the the special case.
   % not (grid has a 1 in either col or row
   %      dist is block
   %      no overlap )
   if ( not( ( Bgrid_spec(1) == 1 || Bgrid_spec(2) == 1) ...
            && strcmp(Bdist_spec(1).dist, 'b') ...
            && strcmp(Bdist_spec(2).dist, 'b') ...
            && isempty(Boverlap) ) )
      % Revert to default.
      A(:,:) = B;

   % Optimized row to column or column to row redistribution
   else
      Alocal = local(A);
      Blocal = local(B);

      % Compute send and receive order
      % Need to shuffle this for best perf.?????
      my_send_order =        circshift(Bproc_list, [0 -(my_rank+1)]);
      my_recv_order = fliplr(circshift(Bproc_list, [0 -(my_rank+1)]));

      pMATLAB.tag_num = pMATLAB.tag_num+1;
      pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

      % Column to row redistribution
      if (Bgrid_spec(1) == 1)
         % Get global ranges of dmats.
         A_Iranges  = global_block_ranges(A,1);
         B_Jranges  = global_block_ranges(B,2);

         % Send out data.
         for dest_rank = my_send_order
            % Get indices.
            i1 = A_Iranges(dest_rank + 1, 2);
            i2 = A_Iranges(dest_rank + 1, 3);
            j1 = B_Jranges(my_rank + 1, 2);
            j2 = B_Jranges(my_rank + 1, 3);
            if (dest_rank == my_rank)
               Alocal(:,j1:j2) = Blocal(i1:i2,:);
            else
               MPI_Send(dest_rank, pMATLAB.tag, comm, Blocal(i1:i2,:));
            end
         end

         % Receive data.
         for recv_rank = my_recv_order
            if (recv_rank ~= my_rank)
               j1 = B_Jranges(recv_rank + 1, 2);
               j2 = B_Jranges(recv_rank + 1, 3);
               Alocal(:,j1:j2) = MPI_Recv(recv_rank, pMATLAB.tag, comm);
            end
         end

      % Row to column redistribution
      else
         % Get global ranges of dmats.
         A_Jranges  = global_block_ranges(A,2);
         B_Iranges  = global_block_ranges(B,1);

         % Send out data.
         for dest_rank = my_send_order
            % Get indices.
            j1 = A_Jranges(dest_rank + 1, 2);
            j2 = A_Jranges(dest_rank + 1, 3);
            i1 = B_Iranges(my_rank + 1, 2);
            i2 = B_Iranges(my_rank + 1, 3);
            if (dest_rank == my_rank)
               Alocal(i1:i2,:) = Blocal(:,j1:j2);
            else
               MPI_Send(dest_rank,pMATLAB.tag,comm,Blocal(:,j1:j2));
            end
         end

         % Receive data.
         for recv_rank = my_recv_order
            if (recv_rank ~= my_rank)
               i1 = B_Iranges(recv_rank + 1, 2);
               i2 = B_Iranges(recv_rank + 1, 3);
               Alocal(i1:i2,:) = MPI_Recv(recv_rank, pMATLAB.tag, comm);
            end
         end
      end

      % Put local data back.
      A = put_local(A, Alocal);

   end
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
