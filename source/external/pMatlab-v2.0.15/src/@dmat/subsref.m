function b = subsref(a, s)
%SUBSREF Subscripted reference. Called for syntax A(S).
%   Should not be called directly.
%   SUBSREF(A, S) Subscripted reference on a distributed array A. 
%       S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts. 
%   !!!WARNING: Does not produce a stand-alone distributed array.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

global pMATLAB;

subs = s(1).subs;
sizeA = size(a);

%
% Array access.
%
if s(1).type == '()' %subscripting type
  % TODO eventually support < cases
  if length(subs) > ndims(a)
    error('@dmat/subsref: Too many dimensions');
  elseif length(subs) < ndims(a)
    error('@dmat/subsref: Too few dimensions');    
  end

  %set submat flag to 0
  submat_flag = 0;
  
  % check to make sure that the indices before the last are
  % within bounds
  f_all = 1;
  for i=1:ndims(a)
    if subs{i} ~= ':'
      f_all = 0;
      if length(subs{i}) ~= 1
          %adjust submat flag
          submat_flag = 1;
      elseif subs{i} > sizeA(i)% && i <= ndims(a)
        error('@dmat/subsref: Index exceeds dmat dimensions');
      end
    end
  end
  
  if ~submat_flag %if reference consist of combinations of : and single numbers, use this code
      % expand the dimensions if needed
      % TODO doesn't follow Matlab semantics if last specified
      % dimension is : and it needs to be expanded, so it has been
      % disabled for the time being
      %
      %excess = subs{length(subs)};
      %for i=length(subs):ndims(a)
      %  subs{i} = mod(excess-1, sizeA) + 1;
      %  excess = floor((excess-1) / sizeA(length(subs))) + 1;
      %end

      %
      % If all subscripts are :, return the entire matrix; otherwise
      % call submatrix.
      %
      if f_all
          b = a;
      else
          %
          % Commonly used variables.
          %
          m = a.map;
          gridA = m.grid;
          distA = m.dist_spec;
          sizeB = sizeA; % initialize size of B to be size of A

          %
          % TODO should use FALLS structure if handling subscript ranges, but
          % currently that is not supported.  Thus using simpler approach.
          %
          for i=1:length(subs)
              if length(subs{i}) == 1
                  if subs{i} ~= ':'
                      if subs{i} < 1 || subs{i} > sizeB(i)
                          error(['@dmat/subsref: The ' num2str(i) 'th subscript ' ...
                              'exceeds size of dmat']);
                      end

                      % figure out distribution
                      switch distA(i).dist
                          case 'b',
                              % block - take the subscript and divide by the block size =
                              % size(a, i) / size(gridA, i)
                              b_size = ceil(size(a, i) / size(gridA, i));
                              idx = floor((subs{i}-1) / b_size) + 1;
                              off = mod(subs{i}-1, b_size) + 1;

                          case 'c',
                              % cyclic - find the remainder of the subscript divided by
                              % the number of processors in that dimension
                              idx = mod(subs{i}-1, size(gridA, i)) + 1;
                              off = floor((subs{i}-1) / size(gridA, i)) + 1;

                          case 'bc',
                              % block cyclic - find out which block this would lie on (block),
                              % and then find out which processor owns this block (cyclic)
                              idx = floor((subs{i}-1) / distA(i).b_size);
                              off = distA(i).b_size * floor(idx / size(gridA, i)) + ...
                                  mod(subs{i}-1, distA(i).b_size) + 1;
                              idx = mod(idx, size(gridA, i)) + 1;

                          otherwise,
                              error(['@dmat/subsref: Unsupported distribution type: ' distA.type]);
                      end

                      %
                      % Set up the subscripts.
                      %
                      s_map.subs{i} = idx;
                      s_data.subs{i} = off;
                      sizeB(i) = 1;
                  else
                      s_map.subs{i} = ':';
                      s_data.subs{i} = ':';
                  end
              else
                  error(['@dmat/subsref: Unsupported subscript: ' subs{i}]);
              end
          end

          %
          % Find the map that would contain these processors and create a
          % dmat using this map.
          %
          s_map.type = '()';
          s_data.type = '()';
          maxGenDim = max([find(sizeB ~= 1) 2]);
          sizeB = sizeB(1:maxGenDim);
          gridB = subsref(gridA, s_map);
          mapB = map(size(gridB), m.dist_spec(1:maxGenDim), reshape(gridB, 1, []));
          b = dmat(sizeB, mapB);

          %
          % If local processor has data that needs to be sent, send it.
          %
          if prod(size(a.local)) > 0 && inmap(mapB, pMATLAB.my_rank)
              b.local = subsref(a.local, s_data);
          end
      end
  else %make a call to submat with a warning
      warning(['dmat/subsref: Fully functional sibscripted reference' ...
          ' is only supported for indices that consist of combinatioons' ...
          ' of : and single number. Otherwise, please restrict operation'...
          ' to the local part of the referenced structure. Stand alone' ...
          ' distirbuted array will not be returned.']);
      b = submat(a,s);
  end
      
  

%
% Structure reference.
%
elseif s(1).type == '.' %subscripting type
    switch subs
        case {'local'}
            b = a.local;
        case{'map'}
            b = a.map;
        otherwise
            error(['@dmat/subsref: ' subs ' cannot be accessed directly or ' ...
                'is not a field of DMAT.']);
    end
end %subscripting type

%
% Recursive call
%
if length(s) > 1
  b = subsref(b, s(2:end));
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
