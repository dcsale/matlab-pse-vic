function p=map(varargin)
%MAP class constructor
%   MAP(GRID_SPEC, DIST_SPEC, PROC_LIST, OVERLAP_SPEC)
%       GRID_SPEC - array of integers specifying how each dimension of a
%           distributed object is broken up.
%           For example is GRID_SPEC = [2 3], the first dimension is broken
%           up between 2 processors and the second dimension is broken up
%           between 3 processors. 
%       DIST_SPEC - array of structures with two possible fields specifying 
%           the distributed array distribution. Each entry in the array has to 
%           have the DIST field defined. The DIST field can have the
%           following values: 
%               b = block
%               c = cyclic
%               bc = block-cyclic
%           Additionally, if DIST == 'bc', the block size 'B_SIZE' must
%           also be defined.
%       PROC_LIST - array of processor ranks specifying on which ranks the
%           object should be distributed. 
%   
%   Returns p:
%   p.DIM       = the number of dimensions of the map (must correspond to the
%       number of dimensions of the the distributed object
%   p.PROC_LIST = the list of processors on which the object should be
%       distributed (varargin{3})
%   p.DIST_SPEC = the distribution description for each dimension
%       (varargin{2})
%   p.GRID      = p.DIM-dimensional array of processors corresponding to how the
%       object should be distributed
%
% Author:   Nadya Travinin

numargs = length(varargin);

%map constructor can be called with or without pMatlab/MatlabMPI
%if the maps are created within the scope of a pMatlab program
%(exist('pMATLAB')), then the processor list is checked against current
%comm scope
if exist('pMATLAB', 'var')==1 
    global pMATLAB;
    Ncpus = pMATLAB.comm_size;
end

%
% Copy constructor.
%
if numargs == 1
  old_map = varargin{1};
  if isstruct(old_map)
    p.dim = old_map.dim;
    p.proc_list = old_map.proc_list;
    p.dist_spec = old_map.dist_spec;
    p.grid = old_map.grid;
    p.overlap = old_map.overlap;  
  elseif isa(old_map, 'map')
    p = old_map;
    return
  else
    error(['MAP CONSTRUCTOR: Using single parameter version requires ' ...
           'a struct with the same fields as a map or a map ' ...
           'itself.']);
  end
  
elseif (numargs == 3) %MAP(GRID_SPEC, DIST_SPEC, PROC_LIST)
    grid_spec = varargin{1};
    dist_spec = varargin{2};
    proc_list = varargin{3};
    
    %ensure that distribution is specified
    if isempty(dist_spec) %default distribution is block
        for d = 1:length(grid_spec)
            dist_spec(d).dist='b';
        end
    elseif length(dist_spec) == 1 
        %if only one distribution is provided, then all dimensions are
        %distributed that way
        if isstr(dist_spec)
            % 'b' for block & 'c' for cyclic distribution
            clear dist_spec;
            for d = 1:length(grid_spec)
                dist_spec(d).dist=varargin{2};
            end
        else
            for d = 2:length(grid_spec)
                dist_spec(d) = dist_spec(1);
            end
        end
    end
    
    dim = length(grid_spec); %dimension of the distributed object
    p.dim = dim;
    p.proc_list = proc_list;
    for i = 1:dim
        %check that distributions defined are consistent with {'b', 'c', 'bc'}
        if (dist_spec(i).dist ~= 'b')&&(dist_spec(i).dist ~= 'c')&&(dist_spec(i).dist ~= 'bc')
            str = strcat('MAP CONSTRUCTOR: ', dist_spec(i).dist, ' is not a valid distribution');
            error(str);
        else   
            %check that block size is defined for block-cyclic distributions
            if (strcmp(dist_spec(i).dist, 'bc'))
                if ~isfield(dist_spec, 'b_size') || isempty(dist_spec(i).b_size)
                    error('MAP CONSTRUCTOR: Block size must be specified for block-cyclic distibution');
                end
            end
        end
    end
    p.dist_spec = dist_spec;
    %create the grid from the processor list
    grid = zeros(grid_spec);
    
   
        %check that the length of the processor list matches the size of
        %the grid
        if (length(proc_list) ~= prod(size(grid)))
            error('MAP CONSTRUCTOR: Processor list does not match the size of the grid');
        else
            grid(:) = proc_list(:);
        end

        %if the maps are created within the scope of a pMatlab program
        %(exist('pMATLAB')), then the processor list is checked against current
        %comm scope
        if exist('pMATLAB', 'var')==1
            %check that the length of the processor list matches the number of
            %processors requested
            if (length(proc_list) > Ncpus)
                warning(['MAP CONSTRUCTOR: Processor list contains more processors ' ...
                    'than number of processors requested.']);
            end
        end
    
    p.grid = grid;    
    p.overlap = []; %no overlap specification
elseif numargs==4   %MAP(GRID_SPEC, DIST_SPEC, PROC_LIST, OVERLAP_SPEC)
    grid_spec = varargin{1};
    dist_spec = varargin{2};
    proc_list = varargin{3};
    overlap_spec = varargin{4};
    
    %ensure that distribution is specified
    if isempty(dist_spec) %default distribution is block
        for d = 1:length(grid_spec)
            dist_spec(d).dist='b';
        end
    elseif length(dist_spec) == 1 
        %if only one distribution is provided, then all dimensions are
        %distributed that way
        if strcmp(dist_spec(1).dist, 'b')
            for d = 2:length(grid_spec)
                dist_spec(d) = dist_spec(1);
            end
        else
            error('MAP CONSTRUCTOR: Overlap is only supported for block distributions.');
        end
    end
    
    dim = length(grid_spec); %dimension of the distributed object
    p.dim = dim;
    p.proc_list = proc_list;
    for i = 1:dim
        %check that distributions defined are consistent with {'b', 'c', 'bc'}
        if ~(strcmp(dist_spec(i).dist, 'b')) 
            error('MAP CONSTRUCTOR: Overlap is only supported for block distributions.');
        end
    end
    
    p.dist_spec = dist_spec;
    %create the grid from the processor list
    grid = zeros(grid_spec);
    %check that the length of the processor list matches the size of
    %the grid
    if (length(proc_list) ~= prod(size(grid)))
        error('MAP CONSTRUCTOR: Processor list does not match the size of the grid');
    else
        grid(:) = proc_list(:);
    end
    p.grid = grid;    
    if (length(overlap_spec) ~= p.dim)
        error('MAP CONSTRUCTOR: Overlap must be specified for all of the dimensions of the map.');
    end
    
    %if the maps are created within the scope of a pMatlab program
    %(exist('pMATLAB')), then the processor list is checked against current
    %comm scope
    if exist('pMATLAB', 'var')==1
        %check that the length of the processor list matches the number of
        %processors requested
        if (Ncpus ~= length(proc_list))
            warning(['MAP CONSTRUCTOR: Processor list does not match the ' ...
                'number of processors requested.']);
        end
    end
    
    p.overlap = overlap_spec;
else
    error('MAP CONSTRUCTOR: Incorrect number of inputs');
end


p=class(p,'map');

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
