function a = subsasgn(a, s, b)
%SUBSASGN Subscripted assignment.
%   A.FIELD = B - allows the fields of a MAP objects to be assigned using
%   the '.' notation (complies with structure behabior)
%
%   This functionality might be deprecated from the final API, to limit
%   control the user has of private members of the MAP object.
%
% Author:  Nadya Travinin
% Edited:  Edmund L. Wong (elwong@ll.mit.edu)

final = length(s) == 1;

if s(1).type=='.'
  switch s(1).subs
   case {'dim'}
    if final
      a.dim = b;
    else
      a.dim = subsasgn(a.dim, s(2:end), b);
    end
   case {'proc_list'}
    if final
      a.proc_list = b;
    else
      a.proc_list = subsasgn(a.proc_list, s(2:end), b);
    end
   case {'dist_spec'}
    if final
      a.dist_spec = b;
    else
      a.dist_spec = subsasgn(a.dist_spec, s(2:end), b);
    end
   case {'grid'}
    if final
      a.grid = b;
    else
      a.grid = subsasgn(a.grid, s(2:end), b);
    end      
   otherwise
    error(strcat(s.subs, ' is not a field of MAP.'));
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