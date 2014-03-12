function [s, varargout] = size(d, dims)
%SIZE Size of the distributed array.
%   SIZE(D) returns the length of each dimension in distributed
%   matrix D.
%
%   SIZE(D, DIMS) returns the length of those dimensions specified
%   in DIMS.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

% if no dimensions are specified, all are used
if nargin == 1
  dims = 1:length(d.size);
end

% If there are no output arguments or 1 output argument,
% return the size as an array
if (nargout == 0 || nargout == 1)
  s(1:length(dims)) = d.size(dims);
else
  s = d.size(dims(1));
  for i = 2:length(dims)
    varargout(i-1) = {d.size(dims(i))};
  end
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
