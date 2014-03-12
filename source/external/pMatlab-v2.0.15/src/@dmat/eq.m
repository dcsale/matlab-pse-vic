function c = eq(a,b)
% == Equal (distributed matrices).
%    A==B compares dimensions, maps and data. 
%    If maps are equal and dimensions agree then the output is a
%    distribituted array with 0 where elements are not equal and 1s where
%    elements are equal (similarly to serial MATLAB).
%    
%    If the maps are not equal, then a 0 is returned regardless of any
%    other properties. Also, if the dimensions and/or sizes are not equal, 
%    then an error gets thrown (analogous to serial MATLAB behavior). 
%
% Author:   Nadya Travinin

%   NOTES:
%    A possible other approach to the case where the maps are not equal, 
%    but the data is, is to return a distributed matrix (distributed along 
%    A's or B's map with 0 and 1 as described above. However, this approach
%    would incur sever communication costs, additionally it is not clear
%    what the distribution of the output matrix should be. Thus for now,
%    this is left as an unimplemented open ended question.

if (a.map==b.map) %compare maps
    if (a.dim==b.dim) %compare dimensions
        if (a.size==b.size) %compare size
            if a.dim==2
                c = ones(a.size(1),a.size(2), a.map);
            elseif a.dim==3
                c = ones(a.size(1), a.size(2), a.size(3), a.map);
            elseif a.dim==4
                c = ones(a.size(1), a.size(2), a.size(3), a.size(4), a.map);
            else
                error('@dmat/eq: Only distributed arrays of up to 4D are supported');
            end
            c.local = (a.local==b.local);
        else %sizes not equal
            error('@dmat/eq:Matrix dimensions must agree');
        end
    else %dimensions not equal
        error('@dmat/eq:Matrix dimensions must agree');
    end
else %maps not equal
    c = logical(0);
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
