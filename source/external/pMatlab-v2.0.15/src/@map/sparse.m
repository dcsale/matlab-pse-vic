function s = sparse(varargin)
%SPARSE Create a sparse distributed matrix.
%    S = SPARSE(I,J,S,M,N,NZMAX,P) generates an M-by-N sparse distributed
%    matrix distributed according map P with space allocated for NZMAX
%    nonzeros (note that NZMAX applies to the overall distributed matrix,
%    not to individual processors.  NZMAX will be distributed as evenly
%    as possible over all processors).
%
%    The rows of [I,J,S] are intended to use to initialize the non-zero
%    values of the matrix.  However, SPARSE currently does not support the
%    use of I, J and S; they are kept to remain consistent with the SPARSE
%    function built into Matlab.
%
%    There are four ways that SPARSE can be called:
%
%    S = SPARSE([],[],[],M,N,NZMAX,P)
%
%    S = SPARSE([],[],[],M,N,P) uses NZMAX = 0.
%  
%    S = SPARSE([],[],[],P) uses M = 0 and N = 0.  This generates the
%    ultimate sparse matrix, an M-by-N all zero matrix.
%  
%    S = SPARSE(M,N,P) abbreviates SPARSE([],[],[],M,N,0,P).  This also
%    generates an M-by-N all zero matrix.
%  
%    The recommended method of creating a sparse distributed matrix is
%    with SPALLOC.
%
%    See also SPALLOC
%
% Author:   Hahn Kim

% Global vars
global pMATLAB;

% SPARSE(M,N,P)
if nargin==3
   ii    = [];          % Row indices for values stored in s
   jj    = [];          % Column indices for values stored in s
   s     = [];          % Non-zero values used to initialize the distributed
                        % sparse matrix
   m     = varargin{1}; % Number of rows
   n     = varargin{2}; % Number of cols
   nzmax = 0;           % Number of non-zero values to allocate space for
   p     = varargin{3}; % Map

% SPARSE(I,J,S,P)
elseif nargin==4
   ii    = varargin{1};
   jj    = varargin{2};
   s     = varargin{3};
   m     = max(ii);
   n     = max(jj);
   nzmax = length(s);
   p     = varargin{4};

% SPARSE(I,J,S,M,N,P)
elseif nargin==6
   ii    = varargin{1};
   jj    = varargin{2};
   s     = varargin{3};
   m     = varargin{4};
   n     = varargin{5};
   nzmax = length(s);
   p     = varargin{6};

% SPARSE(I,J,S,M,N,NZMAX,P)
elseif nargin==7
   ii    = varargin{1};
   jj    = varargin{2};
   s     = varargin{3};
   m     = varargin{4};
   n     = varargin{5};
   nzmax = varargin{6};
   p     = varargin{7};
   
else
   error('@map/sparse: Incorrect number of inputs');
end

% Check if p is of the 'map' class
if ~isa(p, 'map')
   error('@map/spalloc: At least 1 argument must be a map');
end

% For now, sparse only works when ii, jj, and s are [].  Specifying
% values for ii, jj and s will be implemented in the future.
if (length(ii) ~= 0 || length(jj) ~= 0 || length(s) ~= 0)
   error(['@map/sparse: Specifying initial values for sparse '...
          'distributed matrices is not supported, yet.']);
end

% Create the 2D distributed object
s = dmat(m, n, p);

% Figure out local dimensions of dmat

% NOTE: This is recomputing information already computed within
% @dmat/dmat. Is there a cleaner way of getting this information?
falls = get_local_falls(pitfalls(s), p.grid, pMATLAB.my_rank);
local_size = localdims(falls, p.dim);

% What if the user specifies values for ii, jj and s AND a value
% for nzmax, such that on a processor's local portion of the sparse
% dmat, nzmax_local is too small?  Maybe only do this when values
% for i,j and s are specified?  Otherwise, evenly distribute nzmax?
% For now, evenly distribute nzmax across all processors
nzmax_local = floor(nzmax / pMATLAB.comm_size);
nzmax_rem   = nzmax - nzmax_local * pMATLAB.comm_size;
if (pMATLAB.my_rank < nzmax_rem)
   nzmax_local = nzmax_local + 1;
end

% Allocate a sparse matrix for the local portion of the dmat
s.local = spalloc(local_size(1), local_size(2), nzmax_local);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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