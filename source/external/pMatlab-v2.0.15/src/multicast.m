function data = multicast(src, dst, data)
% MULTICAST  Sends data from a single source to multiple destinations.
%    DATA = MULTICAST(SRC, DST, DATA) sends DATA from SRC to DST.
%    DATA will be returned to all destinations as well as the
%    source.  Note that SRC and DST cannot share any common ids;
%    that is, INTERSECT(SRC, DST) = [].
%
%    This method uses O(log2(n)) messages, where n = length(DST).
%
% Author: Edmund Wong (elwong@ll.mit.edu)

global pMATLAB;

% Tag management.
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

% Quick check
if isempty(src) || isempty(dst) || isequal(src, dst)
  return;
end

% Initialize needed variables.
gap = 1;
maxGap = length(dst)/2;
dst = [src reshape(dst, 1, [])];
myIndex = find(dst == pMATLAB.my_rank);

% Communications pattern:
%  1. dst(1)                    => dst(2)
%  2. dst(1), dst(2)            => dst(3), dst(4)
%  ...
%  n. dst(1), ..., dst(2^(n-1)) => dst(2^(n-1)+1), ..., dst(2^n)
%
% Do this until the gap reaches more than halfway across the
% destination list.
while gap <= maxGap
  % send from dst(i) to dst(i+gap)
  if myIndex <= gap
    MPI_Send(dst(myIndex+gap), pMATLAB.tag, pMATLAB.comm, data);
  elseif myIndex <= 2*gap
    data = MPI_Recv(dst(myIndex-gap), pMATLAB.tag, pMATLAB.comm);
  end
  gap = gap*2;
end

% Finish off any remaining destinations.
if myIndex <= length(dst)-gap
  MPI_Send(dst(myIndex + gap), pMATLAB.tag, pMATLAB.comm, data);
elseif myIndex > gap
  data = MPI_Recv(dst(myIndex - gap), pMATLAB.tag, pMATLAB.comm);
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

