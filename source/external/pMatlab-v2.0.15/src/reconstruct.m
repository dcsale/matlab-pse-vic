function mat = reconstruct(pitfalls, grid, temp_mat, mat_size);
%RECONSTRUCT Given collected distributed data in grid layout, reconstructs
%   the original object according to data distribution. If the original
%   distribution was BLOCK, the grid layout is the original layout.
%
%   Also returns reconstructed matrix in the same format as the distributed
%   matrix, i.e. full distributed matrices are returned as full matrices
%   and sparse distributed matrices are returned as sparse matrices.
%
%       MAT_SIZE is the array describing the dimensions of the distributed
%       object.
%
% Author:   Nadya Travinin

%dimension of the distributed object
dim = length(pitfalls);

len = size(grid);
if length(len)<dim
    for i = (length(len)+1):dim
        len(i) = 1;
    end
end

% Determine if the matrix to be reconstructed is a sparse or dense matrix
if (issparse(temp_mat{1,1}))
   % At this point we do not know how many non-zero values the overall
   % sparse distributed matrix contains, so we set it to 0 and let
   % Matlab take care of allocating memory  
   mat = spalloc(mat_size(1), mat_size(2), 0);
else
   % Get the datatype of temp_mat
   datatype = class(temp_mat{1,1});

   % Create a zeros array of the same type as temp_mat
   % Determine Matlab version
%   matlab_version_num = str2num(version('-release'));
%   if (matlab_version_num >= 14)
      mat = zeros(mat_size, datatype);
%   else
%      cmd = ['mat = ' datatype '(zeros(mat_size));'];
%      eval(cmd);
%   end
end

if dim==2
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            local_falls = get_local_falls(pitfalls, grid, grid(i,j));
            global_ind{i,j} = get_global_ind(local_falls, len);
        end
    end
    
    %!!!FOR NOW - assume MAT is of type DOUBLE
    for i = 1:len(1)
        for j = 1:len(2)
          if ~isempty(temp_mat{i, j})
            mat(global_ind{i,j}{1}, global_ind{i,j}{2}) = temp_mat{i,j};
          end
        end
    end    
    
elseif dim==3
    %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                local_falls = get_local_falls(pitfalls, grid, grid(i,j,k));
                global_ind{i,j,k} = get_global_ind(local_falls, len);
            end
        end
    end
    
    %!!!FOR NOW - assume MAT is of type DOUBLE
    for i = 1:len(1)
        for j = 1:len(2)
            for k = 1:len(3)
              if ~isempty(temp_mat{i, j, k})
                mat(global_ind{i,j,k}{1}, global_ind{i,j,k}{2}, global_ind{i,j,k}{3}) = temp_mat{i,j,k};
              end
            end
        end
    end 
elseif dim==4
     %get local indices for each processor in the grid
    for i = 1:len(1)
        for j=1:len(2)
            for k = 1:len(3)
                for m = 1:len(4)
                    local_falls = get_local_falls(pitfalls, grid, grid(i,j,k,m));
                    global_ind{i,j,k,m} = get_global_ind(local_falls, len);
                end
            end
        end
    end
    
    %!!!FOR NOW - assume MAT is of type DOUBLE
    for i = 1:len(1)
        for j = 1:len(2)
            for k = 1:len(3)
                for m = 1:len(4)
                  if ~isempty(temp_mat{i, j, k, m})
                    mat(global_ind{i,j,k,m}{1}, global_ind{i,j,k,m}{2}, ...
                        global_ind{i,j,k,m}{3}, global_ind{i,j,k,m}{4}) = temp_mat{i,j,k,m};
                  end
                end
            end
        end
    end 
else
    error('RECONSTRUCT: Only objects up to 4-D are supported');
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
