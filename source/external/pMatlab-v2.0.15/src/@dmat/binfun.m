function r=binfun(fhandle, p, q)
% Binary operation on matrices P and Q.
%   R = @fhandle(P,Q). Performs the binary operation specified by the 
%   function handle, fhanlde. If both P and Q are distributed, their maps
%   must be the same, otherwise the function throws an error. Additionally,
%   the ditributed arrays is a scalar and the non-distributed isn't, the
%   function throws an error. The addition of a scalar DMAT to a non-scalar
%   non-distributed array would incur significant communication overhead. 
%
% Author:  Nadya Travinin

if isa(p, 'dmat') & ~isa(q, 'dmat')
    %Case 1: First argument (P) is a DMAT.
    %        Second argument (Q) not a DMAT.
    if size(q) == [1 1] 
        %Subcase A: Binary operation with scalar and DMAT.
        p.local = fhandle(p.local,q);
    elseif size(p)==size(q)
        %Subcase B: Dimensions of P and Q are the same.
        my_inds = p.global_ind;
        s = p.size;
        for i = 1:p.dim
             if isa(my_inds{i}, 'char') & (my_inds{i} == ':')
                 my_inds{i} = 1:s(i);
             end
        end    
         p.local = fhandle(p.local, q(my_inds{:}));
    else %Array dimensions do not agree.
        error('dmat/binfun:Array dimensions must agree and the DMAT must be not a scalar.');
    end
    r=p;
elseif ~isa(p, 'dmat') & isa(q, 'dmat')
    %Case 2: First argument (P) is not a DMAT.
    %        Second argument (Q) is a DMAT.
    if size(p) == [1 1]
        %Subcase A: Binary operation with scalar and DMAT.
        q.local = fhandle(p, q.local);
    elseif size(p)==size(q)
        %Subcase B: Dimensions of P and Q are the same.
        my_inds = q.global_ind;
        s = q.size;
        for i = 1:q.dim
            if isa(my_inds{i}, 'char') & (my_inds{i} == ':')
                my_inds{i} = 1:s(i);
            end
        end    
        q.local = fhandle(p(my_inds{:}), q.local);
    else %Array dimensions do not agree
       error('dmat/binfun:Array dimensions must agree and the DMAT must be not a scalar.');
    end
    r=q;
elseif isa(p,'dmat') &  isa(q, 'dmat')
    %Case 3: P and Q are both DMATs.
    % check that the maps are the same
    if p.map ~= q.map
        error('dmat/binfin:Binary operations on distributed arrays are only supported if the maps are the same.');
        %q = remap(q, p.map);
    end
    if size(p)==size(q)
        % add the local parts together
        p.local = fhandle(p.local, q.local);
    else
        error('dmat/binfun: Array dimensions must agree.');
    end
  % copy dmat over
  r = p;
else
    error('dmat/binfun: Unsupported binary operation.');
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
