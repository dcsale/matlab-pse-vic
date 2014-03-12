function OUT = dot(A, B)
% DESCRIPTION:
%   Returns the scalar product of the vectors A and B.
%   A and B must be vectors of the same length.  When A and B are both
%   column vectors, DOT(A,B) is the same as A'*B.
% 
%   For N-D arrays A and B, returns the scalar product
%   along the first non-singleton dimension of A and B. A and B must
%   have the same size.
%
% PARAMETERS:
%   A - first input matrix
%   B - second input matrix
%
% AUTHOR:
%   Daniel Jennings
%   dsj@ll.mit.edu
%
% DATE CREATED:
%   May 31, 2005


global pMATLAB;

%make sure inputs are the same dimension
if length(size(A)) ~= length(size(B))
     error('dmat/dot.m: Inputs A and B must have the same dimensions');
end

%make sure inputs are the same size
if size(A) ~= size(B)
    error('dmat/dot.m: Inputs A and B must be the same size');
end

%check that dimensions are 2-D or 3-D
if (length(size(A)) ~= 2) & (length(size(A)) ~= 3)
    error('dmat/dot.m: Inputs must be either 2-D or 3-D');
end


%*************************************************
%--> adjust all maps so they are along columns <--
%*************************************************
if isa(A, 'dmat')
    %find expected grid spec (want along columns)
    flag = 0;
    if A.dim == 2
        new_gridspec = [1 length(A.map.proc_list)];
        
        %if map and new gridspec do not match up then need remap
        if (size(A.map.grid, 1) ~= new_gridspec(1)) | ...
                (size(A.map.grid, 2) ~= new_gridspec(2))
            flag = 1;
        end
    else
        new_gridspec = [1 length(A.map.proc_list) 1];
        
        %if map and new gridspec do not match up then need remap
        if (size(A.map.grid, 1) ~= new_gridspec(1)) | ...
                (size(A.map.grid, 2) ~= new_gridspec(2)) | ...
                (size(A.map.grid, 3) ~= new_gridspec(3))
            flag = 1;
        end
    end

    %if the grid spec isnt expected then remap
    if flag == 1
        new_proclist = A.map.proc_list;
        
        %TODO is distspec ok?
        new_distspec(1).dist = 'b';
        new_distspec(2).dist = 'b';
        if A.dim == 3
            new_distspec(3).dist = 'b';
        end
        
        newMap = map(new_gridspec, new_distspec, new_proclist);
        A = remap(A, newMap);
    end
end

if isa(B, 'dmat')
    %find expected grid spec (want along columns)
    flag = 0;
    if B.dim == 2
        new_gridspec = [1 length(B.map.proc_list)];
        
        %if map and new gridspec do not match up then need remap
        if (size(B.map.grid, 1) ~= new_gridspec(1)) | ...
                (size(B.map.grid, 2) ~= new_gridspec(2))
            flag = 1;
        end
    else
        new_gridspec = [1 length(B.map.proc_list) 1];
        
        %if map and new gridspec do not match up then need remap
        if (size(B.map.grid, 1) ~= new_gridspec(1)) | ...
                (size(B.map.grid, 2) ~= new_gridspec(2)) | ...
                (size(B.map.grid, 3) ~= new_gridspec(3))
            flag = 1;
        end
    end
    
    %if the grid spec isnt expected then remap
    if flag == 1
        new_proclist = B.map.proc_list;
        
        %TODO is distspec ok?
        new_distspec(1).dist = 'b';
        new_distspec(2).dist = 'b';
        if B.dim == 3
            new_distspec(3).dist = 'b';
        end
        
        newMap = map(new_gridspec, new_distspec, new_proclist);
        B = remap(B, newMap);
    end
end
%*************************************************
%--------> done with all map adjustments <--------
%*************************************************


%check matrix dimensions
if length(size(A)) == 2
%*************************************************
%-----------------> 2-D inputs <------------------
%*************************************************
    if isa(A, 'dmat') & isa(B, 'dmat')    
        %both matrices are distributed
        if A.map ~= B.map
           %make maps equal
           B = remap(B, A.map);
        end

        %maps are equal and already along columns
        OUT = zeros(1, size(A, 2), A.map);
        OUT.local = dot(A.local, B.local);
    elseif isa(A, 'double') & isa(B, 'dmat')            
        %get global indices of B
        b_ind_one = global_ind(B, 1);
        b_ind_two = global_ind(B, 2);
        
        %extract part of A
        atemp = A(b_ind_one, b_ind_two);

        %assign results to OUT
        OUT = zeros(1, size(B, 2), B.map);
        OUT.local = dot(atemp, B.local);
    elseif isa(A, 'dmat') & isa(B, 'double')
        %get global indices of A
        a_ind_one = global_ind(A, 1);
        a_ind_two = global_ind(A, 2);
        
        %extract part of B
        btemp = B(a_ind_one, a_ind_two);

        %assign results to OUT
        OUT = zeros(1, size(A, 2), A.map);
        OUT.local = dot(A.local, btemp);
    else
        %both are non-distributed
        %return normal dot product
        OUT = dot(A, B);
    end
elseif length(size(A)) == 3
%*************************************************
%-----------------> 3-D inputs <------------------
%*************************************************
    if isa(A, 'dmat') & isa(B, 'dmat')    
        %both matrices are distributed
        if A.map ~= B.map
           %make maps equal
           B = remap(B, A.map);
        end

        %maps are equal and already along columns
        OUT = zeros(1, size(A, 2), size(A, 3), A.map);
        OUT.local = dot(A.local, B.local);
    elseif isa(A, 'double') & isa(B, 'dmat')            
        %get global indices of B
        b_ind_one = global_ind(B, 1);
        b_ind_two = global_ind(B, 2);
        b_ind_three = global_ind(B, 3);
        
        %extract part of A
        atemp = A(b_ind_one, b_ind_two, b_ind_three);

        %assign results to OUT
        OUT = zeros(1, size(B, 2), size(B, 3), B.map);
        OUT.local = dot(atemp, B.local);
    elseif isa(A, 'dmat') & isa(B, 'double')
        %get global indices of A
        a_ind_one = global_ind(A, 1);
        a_ind_two = global_ind(A, 2);
        a_ind_three = global_ind(A, 3);
        
        %extract part of B
        btemp = B(a_ind_one, a_ind_two, a_ind_three);

        %assign results to OUT
        OUT = zeros(1, size(A, 2), size(A, 3), A.map);
        OUT.local = dot(A.local, btemp);
    else
        %both are non-distributed
        %return normal dot product
        OUT = dot(A, B);
    end
else
    error('dmat/dot.m: Should never reach this point in the code');
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
