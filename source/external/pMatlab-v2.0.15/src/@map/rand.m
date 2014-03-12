function d = rand(varargin)
%RAND Distributed array of random numbers between 0 and 1 distributed uniformly.
%   NOTE: DIMENSION OF THE DISTRIBUTED ARRAY MUST BE CONSISTENT WITH THE
%       DIMENSION OF THE MAP.
%   Calls the MATLAB RAND function to create each local part of the
%       distributed array. The resulting array will not be the same as a DOUBLE
%       rand array of the same dimensions.
%   RAND(N, P) If N is scalar, an N by N distributed matrix of
%       random numbers mapped according to the map specified by P;
%       if N is a vector, a distributed matrix with dimensions
%       specified by N mapped according to P.
%   RAND(N, P) N by N distributed matrix of random numbers mapped according to the map
%       specified by P.
%   RAND(M, N, P) M by N distributed matrix of random numbers mapped according to the
%       map specified by P.
%   RAND(M, N, Q, P) MxNxQ distributed array of random numbers mapped according to
%       the map specified by P.
%   RAND(M, N, Q, R, P) MxNxQxR distributed array of random numbers mapped according to
%       the map specified by P.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

global pMATLAB;

if nargin<2
  error('map/rand: Number of arguments must be at least 2');
else
  % form dims vector
  dims = varargin{1};
  for i=2:length(varargin)-1
    if length(varargin{i}) > 1
      warning('Input arguments must be scalar.');
    end
    dims = [dims varargin{i}];
  end
  p = varargin{end};
  if ~isa(p, 'map')
    error('map/rand: At least 1 argument must be a map');
  end

  if length(dims) < 5
    d = dmat(dims, p);
  else
    error('map/rand: Incorrect number of inputs');
  end
end

% Allocating memory for the distributed matrix is no longer done
% by @dmat/dmat.  Therefore, we can no longer use the following
% line to get the dimensions of the local portion of the dmat
%local_dims = size(local(d));

% Figure out local dimensions of dmat
% NOTE: This is recomputing information already computed within
% @dmat/dmat. Is there a cleaner way of getting this information?
falls = get_local_falls(pitfalls(d), p.grid, pMATLAB.my_rank);
local_dims = localdims(falls, p.dim);

dim = p.dim;
if dim==2
    g = size(p.grid);
    for j = 1:g(2) %i 1
        for i = 1:g(1)%j 2
            if p.grid(i,j)==pMATLAB.my_rank
                d.local = rand(local_dims);
            else
                rand(local_dims);
            end
        end
    end
elseif dim==3
    [g(1) g(2) g(3)] = size(p.grid);
    for k = 1:g(3)
        for j = 1:g(2)
            for i = 1:g(1)
                if p.grid(i,j,k)==pMATLAB.my_rank
                    d.local = rand(local_dims);
                else
                    rand(local_dims);
                end
            end
        end
    end
elseif dim==4
    [g(1) g(2) g(3) g(4)] = size(p.grid);
    for l = 1:g(4)
        for k = 1:g(3)
            for j = 1:g(2)
                for i = 1:g(1)
                    if p.grid(i,j,k,l)==pMATLAB.my_rank
                        d.local = rand(local_dims);
                    else
                        rand(local_dims);
                    end
                end
            end
        end
    end
else
    error('@MAP/RAND: Only objects up to 4 dimensions are supported.');
end

%d.local = rand(local_dims);

% %>NEW
% dim = p.dim;
% if dim==2
%     for i = 1:n
%         for j = 1:m
%             d(j,i) = rand;
%         end
%     end
% elseif dim==3
%     for k = 1:q
%         for i = 1:n
%             for j = 1:m
%                 d(j,i,k) = rand;
%             end
%         end
%     end
% elseif dim==4
%     for l = 1:r
%         for k = 1:q
%             for i = 1:n
%                 for j=1:m
%                     d(j,i,k,l) = rand;
%                 end
%             end
%         end
%     end
% else
%     error('@MAP/RAND: Only objects up to 4 dimensions are supported.');
% end

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