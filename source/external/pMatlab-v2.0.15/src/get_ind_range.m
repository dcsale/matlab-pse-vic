function ind = get_ind_range(a,s)
%GET_IND_RANGE Creates index ranges from structure S with fields: 
%        type -- string containing '()', '{}', or '.' specifying the
%                subscript type.
%        subs -- Cell array or string containing the actual subscripts. 
%
% Author:   Nadya Travinin

if length(s.subs)==2 %2-D
    ind{1} = s.subs{1};
    ind{2} = s.subs{2};
%     [dims(1) dims(2)] = size(a);
%     if isa(s.subs{1}, 'char') & (s.subs{1}==':') %the first subscipt is a ':'
%         ind{1} = 1:dims(1);
%         if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%             ind{2} = 1:dims(2);
%         else %the second subscript is not ':'
%             ind{2} = s.subs{2};
%         end
%     else  %the first subscript is not ':'
%         ind{1} = s.subs{1};
%         if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%             ind{2} = 1:dims(2);
%         else  %the second subscript is not ':'
%             ind{2} = s.subs{2};
%         end
%     end
elseif length(s.subs)==3 %3-D
    ind{1} = s.subs{1};
    ind{2} = s.subs{2};
    ind{3} = s.subs{3};
%    [dims(1) dims(2) dims(3)] = size(a);
%    if isa(s.subs{1}, 'char') & (s.subs{1}==':') %the first subscipt is a ':'
%        ind{1} = 1:dims(1);
%        if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%            ind{2} = 1:dims(2);
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%            end   %the third subscript is NOT a ':'
%        else %the second subscript is NOT ':'
%            ind{2} = s.subs{2};
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%            end   %the third subscript is NOT a ':'
%        end
%    else    %the first subscipt is NOT a ':'
%        ind{1} = s.subs{1};
%        if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%            ind{2} = 1:dims(2);
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%            end   %the third subscript is NOT a ':'
%        else %the second subscript is NOT ':'
%            ind{2} = s.subs{2};
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%            end   %the third subscript is NOT a ':'
%        end     %the second subscript is NOT ':'
%    end     %the first subscipt is NOT a ':'
elseif length(s.subs)==4 %4-D
    ind{1} = s.subs{1};
    ind{2} = s.subs{2};
    ind{3} = s.subs{3};
    ind{4} = s.subs{4};
%    [dims(1) dims(2) dims(3) dims(4)] = size(a);
%    if isa(s.subs{1}, 'char') & (s.subs{1}==':') %the first subscipt is a ':'
%        ind{1} = 1:dims(1);
%        if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%            ind{2} = 1:dims(2);
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            end   %the third subscript is NOT a ':'
%        else %the second subscript is NOT ':'
%            ind{2} = s.subs{2};
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            end   %the third subscript is NOT a ':'
%        end
%    else    %the first subscipt is NOT a ':'
%        ind{1} = s.subs{1};
%        if isa(s.subs{2}, 'char') & (s.subs{2}==':') %the second subscript is a ':'
%            ind{2} = 1:dims(2);
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            end   %the third subscript is NOT a ':'
%        else %the second subscript is NOT ':'
%            ind{2} = s.subs{2};
%            if isa(s.subs{3}, 'char') & (s.subs{3}==':') %the third subscript is a ':'
%                ind{3} = 1:dims(3);
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            else    %the third subscript is NOT a ':'
%                ind{3} = s.subs{3};
%                if isa(s.subs{4}, 'char') & (s.subs{4}==':') %the fourth subscript is a ':'
%                    ind{4} = 1:dims(4);
%                else %the fourth subscript is NOT ':'
%                    ind{4} = s.subs{4};
%                end %the fourth subscript is NOT ':'
%            end   %the third subscript is NOT a ':'
%        end     %the second subscript is NOT ':'
%    end     %the first subscipt is NOT a ':'
else
    error('GET_IND_RANGE: Only up to 4 dimensional objects are supported.')
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