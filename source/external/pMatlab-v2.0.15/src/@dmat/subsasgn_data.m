function [data, a_local_ind] = subsasgn_data(a, b, falls_index, fi)
%SUBSASGN_DATA Helper function for distributed array subsasgn. 
%   SUBSASGN_DATA(A, B, FALLS_INDEX, FI) Computes which data to send from  
%       distributed array B to distributed array A based on falls
%       intersection FI.
%
% Author:   Nadya Travinin

%dimension of the distributed object
dim = length(fi);

if dim==2
    num_data=0;
    for f1 = 1:length(fi{1}{falls_index})
        for f2 = 1:length(fi{2}{falls_index})
            clear fi_temp;
            fi_temp(1) = fi{1}{falls_index}(f1);
            fi_temp(2) = fi{2}{falls_index}(f2);
            num_data=num_data+1;
            g_ind = get_global_ind(fi_temp);
            b_local_ind = get_local_ind(b.global_ind, g_ind);
            data{num_data} = b.local(b_local_ind{1}, b_local_ind{2});
            a_local_ind{num_data} = get_local_ind(a.global_ind, g_ind);
        end
    end
elseif dim==3
    num_data=0;
    for f1 = 1:length(fi{1}{falls_index})
        for f2 = 1:length(fi{2}{falls_index})
            for f3 = 1:length(fi{3}{falls_index})
                clear fi_temp;
                fi_temp(1) = fi{1}{falls_index}(f1);
                fi_temp(2) = fi{2}{falls_index}(f2);
                fi_temp(3) = fi{3}{falls_index}(f3);
                num_data=num_data+1;
                g_ind = get_global_ind(fi_temp);
                b_local_ind = get_local_ind(b.global_ind, g_ind);
                data{num_data} = b.local(b_local_ind{1}, b_local_ind{2}, b_local_ind{3});
                a_local_ind{num_data} = get_local_ind(a.global_ind, g_ind);
               
            end %f3
        end %f2
    end %f1
elseif dim==4
    num_data=0;
    for f1 = 1:length(fi{1}{falls_index})
        for f2 = 1:length(fi{2}{falls_index})
            for f3 = 1:length(fi{3}{falls_index})
                for f4 = 1:length(fi{4}{falls_index})
                    clear fi_temp;
                    fi_temp(1) = fi{1}{falls_index}(f1);
                    fi_temp(2) = fi{2}{falls_index}(f2);
                    fi_temp(3) = fi{3}{falls_index}(f3);
                    fi_temp(4) = fi{4}{falls_index}(f4);
                    num_data=num_data+1;
                    g_ind = get_global_ind(fi_temp);
                    b_local_ind = get_local_ind(b.global_ind, g_ind);
                    data{num_data} = b.local(b_local_ind{1}, b_local_ind{2}, b_local_ind{3}, b_local_ind{4});
                    a_local_ind{num_data} = get_local_ind(a.global_ind, g_ind);
                   
                end %f4
            end %f3
        end %f2
    end %f1
else
    error('DMAT/SUBSASGN_DATA: Only objects up to four (4) dimensions are supported.');
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
