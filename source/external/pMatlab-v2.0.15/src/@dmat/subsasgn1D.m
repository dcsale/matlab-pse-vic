function a = subsasgn1D(a,s,b)
%SUBSASGN1D One dimensional subsasgn.
%   S is of the following form (:), independent of the dimension of the
%   distributed object dimension.
%
% Author:   Nadya Travinin

global pMATLAB;

% Instead of creating a copy of a, write directly to the memory
% allocated for a in the caller's workspace
assignin('caller', inputname(1), []);

if isa(b, 'double')
    if s.subs{1}==':'
        if (length(size(b))==2) & (size(b)==[1 1]) %b is a scalar
            %if (size(b)==[1 1])
            if inmap(a.map, pMATLAB.my_rank)
                %assigment to a scalar
                a.local(:) = b;
            end
        else
            %b is a regular non distributed matrix
            %check that dimensions are the same and redistribute
            %according to a's map
            if (size(b) == a.size) %dimensions are the same
                if a.dim == 2 %2-D
                    a.local(:,:) = b(a.global_ind{1}, a.global_ind{2});
                elseif a.dim == 3 %3-D
                    a.local(:,:,:) = b(a.global_ind{1}, a.global_ind{2}, a.global_ind{3});
                elseif a.dim == 4 %4-D
                    a.local(:,:,:,:) = b(a.global_ind{1}, a.global_ind{2}, a.global_ind{3}, a.global_ind{4});
                else %dimension > 4
                    error('DMAT/SUBSASGN1D: Only up to 4 dimensional objects supported.');
                end %distributed object dimension
            end %dimensions match
        end
    else
        error('unsupported indexing');
    end
elseif isa(b, 'dmat') %RHS is a DMAT
    if isa(s.subs{1}, 'char')  %subscript is a CHAR
        if s.subs{1} == ':' %subscript is a ':'
            if a.dim == 2 %2-D
                s.subs{2} = ':';
                a = subsasgn2D(a,s,b);
            elseif a.dim == 3 %3-D
                s.subs{2} = ':';
                s.subs{3} = ':';
                a = subsasgn3D(a,s,b);
            elseif a.dim == 4 %4-D
                s.subs{2} = ':';
                s.subs{3} = ':';
                s.subs{4} = ':';
                a = subsasgn4D(a,s,b);
            else %>4-D
                error('DMAT/SUBSASGN1D: Only up to 4 dimensional objects supported.');
            end  %>4-D 
        else  %subscript is not ':'
            error('unsupported indexing');
        end  %subscript is not ':'
    else  %subscript is NOT a CHAR
        error('unsupported indexing');
    end  %subscript is NOT a CHAR
else %RHS is a not a DMAT or a DOUBLE
    error('DMAT/SUBSASGN1D: RHS must be a DOUBLE or DMAT.');
end  %RHS is a not a DMAT or a DOUBLE

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
