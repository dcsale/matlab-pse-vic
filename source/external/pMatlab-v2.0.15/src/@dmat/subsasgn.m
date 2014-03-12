function a = subsasgn(a, s, b)
%SUBSASGN Subscripted assignment to a distributed object. Called for syntax A(I) = B.
%   Should not be called directly.
%   SUBSASGN(A, S, B) Subscripted assignment of B (right hand side) to A
%       (left hand side). A is a DMAT (distributed array). B can be a DMAT or a
%       DOUBLE. S is a structure array with the fields:
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts. 
%
% Author:   Nadya Travinin

% Instead of creating a copy of a, write directly to the memory
% allocated for a in the caller's workspace

%  Following line doesn't work in Octave.
%assignin('caller', inputname(1), []);

if length(s)==1 %subscripting level 
    if s.type == '()' %subscripting type -> parenthesis
        if length(s.subs)==1  %1-D subscripted assignment
            a = subsasgn1D(a, s, b);
        elseif length(s.subs) == a.dim %number of dimensions of indices must match the number of dimensions of the matrix
            if length(s.subs)==2  %2-D subscripted assignment
                a = subsasgn2D(a,s,b);
            elseif length(s.subs)==3 %3-D subscripted assignment
                a = subsasgn3D(a,s,b);
            elseif length(s.subs)==4 %4-D subscripted assignment
                a = subsasgn4D(a,s,b);
            else
                error('DMAT/SUBSASGN: Only objects up to four (4) dimensions are supported.');
            end %N-D subscripted assignment
        else
            error('DMAT/SUBSASGN: The number of index dimensions must match the number of dimensions of the distributed array.');
        end
    elseif s.type == '.' %subscripting type - members of dmat structure
        if s.subs=='local'
            a.local=b;
        end
    end
elseif length(s) == 2  %subscripting level
    if s(1).type == '.' %subscripting type - members of dmat structure
        if s(1).subs == 'local'
            a.local(s(2).subs{:}) = b;
        else
            error('DMAT/SUBSASGN: Incorrect subscripting.');
        end
    else %subscripting type
        error('DMAT/SUBSASGN: Incorrect subscripting.');
    end %subscripting type
else %subscripting level > 2
    error('DMAT/SUBSASGN: Incorrect subscripting.');
end %subscripting level

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
