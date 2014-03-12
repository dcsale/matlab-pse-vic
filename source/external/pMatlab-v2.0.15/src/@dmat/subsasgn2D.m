function a = subsasgn2D(a,s,b)
%SUBSASGN2D Two dimensional subsasgn.
%   S is of the following form (i:j, k:l). Distributed object's dimension
%   is 2. B is either a DMAT or a DOUBLE.
%
% Author:   Nadya Travinin

global pMATLAB;

% Instead of creating a copy of a, write directly to the memory
% allocated for a in the caller's workspace
assignin('caller', inputname(1), []);

if isa(b, 'double') %RHS is a double
    if isa(s.subs{1}, 'char') & isa(s.subs{2}, 'char')
        if (s.subs{1} == ':') & (s.subs{2} == ':') %A(:,:) = B 
            if (size(b) == a.size) %dimensions are the same
                a.local(:,:) = b(a.global_ind{1}, a.global_ind{2});
            elseif size(b)==[1 1]  %b is a scalar
                a.local(:,:) = b;
            else %dimensions do not match
                error('DMAT/SUBSASGN2D:  Subscripted assignment dimension mismatch.');
            end
        end
    else %A(i:j, k:l) = B
        ind = get_ind_range(a,s);
        local_ind = get_local_ind(a.global_ind, ind);
            
        if size(b) == size(a.local(local_ind{1}, local_ind{2}))
            a.local(local_ind{1}, local_ind{2}) = b;
        end

    end %A(i:j, k:l) = B
elseif isa(b, 'dmat') %assignment from a distributed matrix
    %communication might be necessary 
    if (s.subs{1} == ':') & (s.subs{2} == ':') %A(:,:) = B
        %check that dimensions match
        if a.size ~= b.size
            error('DMAT/SUBSASGN2D: Subscripted assignment dimension mismatch.');
        end
        
        %check if maps are the same
        if a.map==b.map %maps are the same - no communication needed
            a.local(:,:) = b.local(:,:);
        else %maps not the same - redistribution
            
            %compute falls intersections
            if (inmap(a.map, pMATLAB.my_rank)) | (inmap(b.map, pMATLAB.my_rank)) 
                %the local processor is either in a's or b's map
                
                %if local processor belongs to b's map, get
                %a's local falls and compute intersections
                if inmap(b.map, pMATLAB.my_rank) %belongs to b's map
                    for i=1:length(a.map.proc_list)
                        a_falls = get_local_falls(a.pitfalls, a.map.grid, a.map.proc_list(i));
                        %falls intersection on b's procs
                        b_row_fi{i} = falls_intersection(b.falls(1), a_falls(1));    
                        b_col_fi{i} = falls_intersection(b.falls(2), a_falls(2));   
                    end       
                    
                end %belongs to b's map
                
                %if local processor belongs to a's map, get
                %b's local falls and compute falls
                %intersections
                if inmap(a.map, pMATLAB.my_rank)
                    for i = 1:length(b.map.proc_list)
                        b_falls = get_local_falls(b.pitfalls, b.map.grid, b.map.proc_list(i));
                        %falls intersection on a's procs
                        a_row_fi{i} = falls_intersection(b_falls(1), a.falls(1));
                        a_col_fi{i} = falls_intersection(b_falls(2), a.falls(2));
                    end
                end %belongs to a's map   
                
            end %the local processor is either in a's or b's map, otherwise should just fall through
            
              %determine which data to send and send the data
        for p1=1:length(b.map.proc_list) %iterate through b's processor list
            for p2=1:length(a.map.proc_list) %iterate through a's processor list
                
                %increment tag
                pMATLAB.tag_num = pMATLAB.tag_num+1;
                pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));
                
                if (inmap(a.map, pMATLAB.my_rank)) | (inmap(b.map, pMATLAB.my_rank)) 
                    %the local processor is either in a's or b's map
                    if b.map.proc_list(p1) ~= a.map.proc_list(p2) %comm is needed
                        if pMATLAB.my_rank==b.map.proc_list(p1) %my rank is current B rank
                            %redistribute data
                            if ~isempty(b_row_fi{p2}) & ~isempty(b_col_fi{p2}) %both intersections not empty
                                %[data, scratch] = subsasgn_data(a, b, p2, b_row_fi, b_col_fi) ;
                                b_fi{1} = b_row_fi;
                                b_fi{2} = b_col_fi;
                                [data, scratch] = subsasgn_data(a, b, p2, b_fi);
                                MPI_Send(a.map.proc_list(p2), pMATLAB.tag, pMATLAB.comm, data);
                            end %both intersections not empty
                        elseif pMATLAB.my_rank==a.map.proc_list(p2) %my_rank is current A rank
                            if ~isempty(a_row_fi{p1}) & ~isempty(a_col_fi{p1}) %both intersections not empty
                                %[scratch, a_local_ind] = subsasgn_data(a, b, p1, a_row_fi, a_col_fi) ;
                                a_fi{1} = a_row_fi;
                                a_fi{2} = a_col_fi;
                                [scratch, a_local_ind] = subsasgn_data(a, b, p1, a_fi) ;
                                data = MPI_Recv(b.map.proc_list(p1), pMATLAB.tag, pMATLAB.comm);
                                for r = 1: length(data)
                                    ind = a_local_ind{r};
                                    a.local(ind{1}, ind{2}) = data{r};
                                end
                            end %both intersections not empty
                        end   %my_rank is current A rank
                    elseif b.map.proc_list(p1) == a.map.proc_list(p2) %no comm 
                        if pMATLAB.my_rank==a.map.proc_list(p2)
                            if ~isempty(b_row_fi{p2}) & ~isempty(b_col_fi{p2}) %both intersections not empty
                                %[data, a_local_ind] = subsasgn_data(a, b, p2, b_row_fi, b_col_fi) ;
                                b_fi{1} = b_row_fi;
                                b_fi{2} = b_col_fi;
                                [data, a_local_ind] = subsasgn_data(a, b, p2, b_fi);
                                for r = 1: length(data)
                                    ind = a_local_ind{r};
                                    a.local(ind{1}, ind{2}) = data{r};
                                end
                            end  %both intersections not empty
                        end
                    end %no comm
                end %the local processor is either in a's or b's map, otherwise should just fall through      
            end %iterate thorugh a's processor list
        end %iterate through b's processor list
            
        end %maps not the same - redistribution
    else  %A(i:j, k:l) = B
        ind{1} = s.subs{1};
        ind{2} = s.subs{2};
        a_local_ind = get_local_ind(a.global_ind, ind);
        if size(a.local(a_local_ind{:}))==size(b.local) 
            a.local(a_local_ind{:}) = b.local;
        else
            error('DMAT/SUBSASGN2D: Subscripted assignment dimensions mismatch.')
        end
    end  %A(i:j, k:l) = B 
else %RHS is not a DMAT or a DOUBLE
     error('DMAT/SUBSASGN2D: RHS must be a DOUBLE or DMAT.');
end  %b is a dmat

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
