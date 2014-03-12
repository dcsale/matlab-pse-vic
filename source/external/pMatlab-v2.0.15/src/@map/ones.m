function o=ones(varargin)
%ONES Ones distributed array.
%   NOTE: DIMENSION OF THE DISTRIBUTED ARRAY MUST BE CONSISTENT WITH THE
%       DIMENSION OF THE MAP.
%   ONES(N, P) If N is scalar, an N by N distributed matrix of
%       ones mapped according to the map specified by P; if N is a
%       vector, a distributed matrix with dimensions specified by N
%       mapped according to P.
%   ONES(M, N, P) M by N distributed matrix of ones mapped according to the
%       map specified by P.
%   ONES(M, N, Q, P) MxNxQ distributed array of ones mapped according to
%       the map specified by P.
%   ONES(M, N, Q, R, P) MxNxQxR distributed array of ones mapped according to
%       the map specified by P.
%   ONES(M, N, ..., P, TYPE) MxNx... distributed array of ones of datatype
%       TYPE mapped according to the map specified by P.
%
%       Example:
%          Create a 100x10 dmat of 8-bit signed integers
%          p = map([1 Ncpus], {}, 0:Ncpus-1);
%          x = ones(100, 10, p, 'int8');
%
% Author:  Edmund L. Wong (elwong@ll.mit.edu)

% Global vars - required to compute falls information
global pMATLAB;

if nargin<2
  error('map/ones: Number of arguments must be at least 2');
else
  % form dims vector
  dims = varargin{1};

  % Check if last argument is a string
  if (strcmp(class(varargin{end}), 'char'))
    datatype = varargin{end};

    % Ignore map and database arguments
    ndims = nargin - 2;

    % Get the map argument
    p = varargin{end - 1};
  else
    datatype = 'double';
  
    % Ignore map argument
    ndims = nargin - 1;

    % Get the map argument
    p = varargin{end};
  end

  for i=2:ndims
    if length(varargin{i}) > 1        
      warning('Input arguments must be scalar.');
    end
    dims = [dims varargin{i}];
  end

  if ~isa(p, 'map')
    error('map/ones: At least 1 argument must be a map');
  end

  if length(dims) < 5
    d = dmat(dims, p);
  else
    error('map/ones: Incorrect number of inputs');
  end
end

% Figure out local dimensions of dmat
% NOTE: This is recomputing information already computed within
% @dmat/dmat. Is there a cleaner way of getting this information?
falls = get_local_falls(pitfalls(d), p.grid, pMATLAB.my_rank);
local_dims = localdims(falls, p.dim);

% Allocating memory for the distributed matrix is no longer done
% by @dmat/dmat.
%d(:)=1;

% Allocate a ones matrix for the local portion of the dmat
% Determine Matlab version
matlab_version_num = str2num(version('-release'));
if (matlab_version_num >= 14)
   d.local = ones(local_dims, datatype);
else
   cmd = ['d.local = ' datatype '(ones(local_dims));'];
   eval(cmd);
end

o = d;

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