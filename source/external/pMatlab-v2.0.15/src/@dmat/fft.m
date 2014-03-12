function x = fft(x, varargin)
%FFT Discrete Fourier transform on a distributed matrix.
%   FFT(X) is the discrete Fourier transform (DFT) of matrix X. The FFT
%   operation is applied to each column. If the matrix X is row
%   distributed, throws a warning and remaps the distributed array. Calls
%   the MATLAB fft on the local part.
%
%   FFT(X,[],DIM) or FFT(X,N,DIM) applies the FFT operation across the
%   dimension DIM. If the matrix is distributed along different dimension
%   than DIM, throws a warning and remaps the array along appropriate
%   dimension. Calls the MATLAB fft on the local part.
%
%   NOTE: If the fft remaps X, returned X has a new map (different from the
%   orginal map passed in).
%
% Author:   Nadya Travinin

%Updates to @dmat/fft
%Changes: 
%   (1) Supports the fft operation on 3-D numerical arrays
%   (2) If the numerical array passed into the fft is not distributed
%   accordingly, remaps the array warning the user that remapping is taking
%   place
%       Issue to consider: Should the array be mapped to the original
%       mapping after the fft is performed? 
%       Plus: The output is mapped the same as the input
%       Minus: Significant overhead due to communication required during
%       remapping

if x.dim == 2  %distributed array is a matrix
    g = grid(x);
    grid_dims = size(g);
    if nargin==1 %FFT(X): default calling convention - fft along the columns
        if grid_dims(1) == 1 %the matrix is broken up along the columns
            x.local = fft(x.local);
        else %the matrix is broken up along the rows, redistribute
            warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along columns.')
            %>>REMAPPING CODE
            grid_spec = [1 grid_dims(1)*grid_dims(2)];
            old_map = x.map;
            dist_spec = old_map.dist_spec;
            proc_list = g(:)';
            new_map = map(grid_spec, dist_spec, proc_list);
            x = remap(x,new_map);
            %>>REMAPPING CODE
            x.local = fft(x.local);
        end
    elseif nargin == 3 %FFT(X,[],DIM)
        if varargin{2}==2 %fft along the rows
            if grid_dims(2)==1
                x.local = fft(x.local, varargin{1}, varargin{2});
            else
                warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along rows.')
                %>>REMAPPING CODE
                grid_spec = [grid_dims(1)*grid_dims(2) 1];
                old_map = x.map;
                dist_spec = old_map.dist_spec;
                proc_list = g(:)';
                new_map = map(grid_spec, dist_spec, proc_list);
                x = remap(x,new_map);
                %>>REMAPPING CODE
                x.local = fft(x.local, varargin{1}, varargin{2});
            end
        elseif varargin{2}==1 %fft applied to each column
            if grid_dims(1) == 1
                x.local = fft(x.local, varargin{1}, varargin{2});
            else
                warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along columns.')
                %>>REMAPPING CODE
                grid_spec = [1 grid_dims(1)*grid_dims(2)];
                old_map = x.map;
                dist_spec = old_map.dist_spec;
                proc_list = g(:)';
                new_map = map(grid_spec, dist_spec, proc_list);
                x = remap(x,new_map);
                %>>REMAPPING CODE
                x.local = fft(x.local, varargin{1}, varargin{2});
            end
        end
    else
        error('@dmat/fft: Distributed fft must be called with either 1 or 3 arguments.');
    end
    
elseif x.dim==3 %distributed array is 3-D
    
    g = grid(x);
    grid_dims = size(g);
    if nargin==1 %FFT(X): default calling convention - fft along the columns
        disp('3D FFT: fft(x)');
        if grid_dims(1) == 1 %the matrix is broken up along the columns & the third dimension, but not rows
            x.local = fft(x.local);
        else %the matrix is broken up along the rows, redistribute
            warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along 3-rd dimension.')
            %>>REMAPPING CODE
            if length(grid_dims)<3 %fill in singleton dimensions
                grid_dims(3)=1;
            end
            grid_spec = [1  1 grid_dims(1)*grid_dims(2)*grid_dims(3)]; %map along the third dimension
            %!!!NEED TO CHECK THAT THIS IS IN FACT THE
            %MOST EFFICIENT MAPPING
            old_map = x.map;
            dist_spec = old_map.dist_spec;
            proc_list = g(:)';
            new_map = map(grid_spec, dist_spec, proc_list);
            x = remap(x,new_map);
            %>>REMAPPING CODE
            x.local = fft(x.local);
        end
    elseif nargin == 3 %FFT(X,[],DIM)
        if varargin{2}==2 %fft along the rows
            disp('3D FFT: fft(x, [], 2)');
            if grid_dims(2) == 1 %the matrix is broken up along the rows & the third dimension, but not columns
                x.local = fft(x.local,varargin{1}, varargin{2});
            else %the matrix is broken up along the columns, redistribute
                warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along 3-rd dimension.')
                %>>REMAPPING CODE
                if length(grid_dims)<3 %fill in singleton dimensions
                    grid_dims(3)=1;
                end
                grid_spec = [1  1 grid_dims(1)*grid_dims(2)*grid_dims(3)]; %map along the third dimension
                %!!!NEED TO CHECK THAT THIS IS IN FACT THE
                %MOST EFFICIENT MAPPING
                old_map = x.map;
                dist_spec = old_map.dist_spec;
                proc_list = g(:)';
                new_map = map(grid_spec, dist_spec, proc_list);
                x = remap(x,new_map);
                %>>REMAPPING CODE
                x.local = fft(x.local,varargin{1}, varargin{2});
            end
        elseif varargin{2}==1 %fft applied to each column 
            disp('3D FFT: fft(x, [], 1)');
            if grid_dims(1) == 1 %the matrix is broken up along the columns & the third dimension, but not rows
                x.local = fft(x.local,varargin{1}, varargin{2});
            else %the matrix is broken up along the rows, redistribute
                warning('@dmat/fft: The matrix is not mapped along the appropriate dimension, remapping along 3-rd dimension.')
                %>>REMAPPING CODE
                if length(grid_dims)<3 %fill in singleton dimensions
                    grid_dims(3)=1;
                end
                grid_spec = [1  1 grid_dims(1)*grid_dims(2)*grid_dims(3)]; %map along the third dimension
                %!!!NEED TO CHECK THAT THIS IS IN FACT THE
                %MOST EFFICIENT MAPPING
                old_map = x.map;
                dist_spec = old_map.dist_spec;
                proc_list = g(:)';
                new_map = map(grid_spec, dist_spec, proc_list);
                x = remap(x,new_map);
                %>>REMAPPING CODE
                x.local = fft(x.local,varargin{1}, varargin{2});
            end
        end
    else
        error('@dmat/fft: Distributed fft must be called with either 1 or 3 arguments.');
    end
    
else
    error('@dmat/fft: FFT can only be applied to matrices or 3-D arrays.');
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
