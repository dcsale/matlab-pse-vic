function r = get_local_proc(pitfalls, grid, ind)
%GET_LOCAL_PROC Returns the rank of the processor that contains index IND.
%   GET_LOCAL_IND(PITFALLS, GRID, IND) Given the pitfalls structure, the
%   grid and index pair, computes the rank of the processor where the index
%   pair is local. IND is an array of two elements, GRID is the processor grid.
%
% Author:   Nadya Travinin

%get dimensions of the grid
g_dims = size(grid);

% 2D processor grid
if (length(g_dims) == 2)
   %get global indices for each processor in the grid
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         local_falls = get_local_falls(pitfalls, grid, grid(i,j));
         global_ind{i,j} = get_global_ind(local_falls);
      end
   end

   %search each processor's global indices for the requested indices
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         if (  ~isempty(find(global_ind{i,j}{1}==ind(1))) ...
             & ~isempty(find(global_ind{i,j}{2}==ind(2))) )
            r = grid(i,j);
         end
      end
   end 

% 3D processor grid
elseif (length(g_dims) == 3)
   %get global indices for each processor in the grid
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         for k = 1:g_dims(3)
            local_falls = get_local_falls(pitfalls, grid, grid(i,j,k));
            global_ind{i,j,k} = get_global_ind(local_falls);
         end
      end
   end

   %search each processor's global indices for the requested indices
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         for k = 1:g_dims(3)
            if (  ~isempty(find(global_ind{i,j,k}{1}==ind(1))) ...
                & ~isempty(find(global_ind{i,j,k}{2}==ind(2))) ...
                & ~isempty(find(global_ind{i,j,k}{3}==ind(3))) )
               r = grid(i,j,k);
            end
         end
      end
   end 

% 4D processor grid
elseif (length(g_dims) == 4)
   %get global indices for each processor in the grid
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         for k = 1:g_dims(3)
            for l = 1:g_dims(4)
               local_falls = get_local_falls(pitfalls, grid, grid(i,j,k,l));
               global_ind{i,j,k,l} = get_global_ind(local_falls);
            end
         end
      end
   end

   %search each processor's global indices for the requested indices
   for i = 1:g_dims(1)
      for j = 1:g_dims(2)
         for k = 1:g_dims(3)
            for l = 1:g_dims(4)
               if (  ~isempty(find(global_ind{i,j,k,l}{1}==ind(1))) ...
                   & ~isempty(find(global_ind{i,j,k,l}{2}==ind(2))) ...
                   & ~isempty(find(global_ind{i,j,k,l}{3}==ind(3))) ...
                   & ~isempty(find(global_ind{i,j,k,l}{4}==ind(4))) )
                  r = grid(i,j,k,l);
               end
            end
         end
      end
   end 

else
   error(['GET_LOCAL_PROC: does not support ' num2str(g_dims) ' dimensions.']);
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