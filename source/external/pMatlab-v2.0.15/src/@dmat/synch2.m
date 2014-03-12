function d = synch2(d)
% SYNCH Syncronize the data in the distribute matrix. 
%   SYNCH(D) No-op if there is no overlap. If overlap is present, the owner
%       processor of the overlaping data sends its data to the processor that
%       owns the copy of the overlapping data. The owner is the processor with
%       the higher index in the grid in the corresponding dimension. For
%       example, if the overlap is in the second dimension the owner is the
%       processor in the column of the grid with the higher index.
%
% Author:   Nadya Travinin
%
% Dec 6, 2005 - added functionality to allow for 2-way shifting of data to
% support 2 way over-lap.  Julie Mullen 

global pMATLAB;

p = d.map;
%find current processor index in the grid
proc_grid = p.grid;
my_grid_ind = n_dim_find(proc_grid, pMATLAB.my_rank);
grid_dims = size(proc_grid);

%initialize the buffers for the data to be passed
data_m1 = cell(1,d.dim);
data_p1 = cell(1,d.dim);
for s = 1:d.dim
    data_m1{s} = ':';
    data_p1{s} = ':';
end

for i = 1:d.dim %syncronize each dimension
%increment tags
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

if inmap(p, pMATLAB.my_rank)
    if ~isempty(p.overlap) %overlap description defined
       if p.overlap(i)>0 %overlap is greater than 0
           
            subs = cell(1,d.dim);
            for s = 1:d.dim
                subs{s} = ':';
            end
            %buffer the overlap data                
             if my_grid_ind(i) ~= 1
                subs{i} = 1:p.overlap(i);
                data_m1{i} = d.local(subs{:});
             end % buffer of data going to my_rank-1
             
             if my_grid_ind(i) ~= grid_dims(i)
                 local_dims = size(d.local);
                 subs{i} = (local_dims(i)-p.overlap(i)+1):local_dims(i); 
                 data_p1{i} = d.local(subs{:});
             end % buffer of data going to my_rank+1
             
             ind_adjuster = zeros(1,d.dim);
             ind_adjuster(i) = -1;
             recv_proc_ind = my_grid_ind+ind_adjuster;
             
                     
             if my_grid_ind(i) ~= 1  
                 if d.dim==2
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2)), pMATLAB.tag, pMATLAB.comm, data_m1{i})
                 elseif d.dim==3
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2), recv_proc_ind(3)), pMATLAB.tag, pMATLAB.comm, data_m1{i});
                 elseif d.dim==4
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2), recv_proc_ind(3), recv_proc_ind(4)), pMATLAB.tag, pMATLAB.comm, data_m1{i});
                 end
             end
                
             ind_adjuster(i) = 1;
             send_proc_ind = my_grid_ind+ind_adjuster;
             
             if my_grid_ind(i) ~= grid_dims(i)
                 if d.dim==2
                     data_m1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2)), pMATLAB.tag, pMATLAB.comm)
                 elseif d.dim==3
                     data_m1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2), send_proc_ind(3)), pMATLAB.tag, pMATLAB.comm);
                 elseif d.dim==4
                     data_m1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2), send_proc_ind(3), send_proc_ind(4)), pMATLAB.tag, pMATLAB.comm)
                 end
             end

           end %overlap is greater than 0
        end %overlap description defined  
    end % inmap
                
%now send data from left->right
% INCREMENT TAGS FOR ALL PROCESSORS
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

if inmap(p, pMATLAB.my_rank)
     if ~isempty(p.overlap) %overlap description defined
         if p.overlap(i)>0 %overlap is greater than 0
                
             ind_adjuster = zeros(1,d.dim);
             %ind_adjuster is 1, sending to my_rank+1
             ind_adjuster(i) = 1;
             recv_proc_ind = my_grid_ind+ind_adjuster;
                
             if my_grid_ind(i) ~= grid_dims(i)                
                 if d.dim==2
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2)), pMATLAB.tag, pMATLAB.comm, data_p1{i});
                 elseif d.dim==3
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2), recv_proc_ind(3)), pMATLAB.tag, pMATLAB.comm, data_p1{i});
                 elseif d.dim==4
                     MPI_Send(proc_grid(recv_proc_ind(1), recv_proc_ind(2), recv_proc_ind(3), recv_proc_ind(4)), pMATLAB.tag, pMATLAB.comm, data_p1{i});
                 end
             end
                
             ind_adjuster(i) = -1;
             send_proc_ind = my_grid_ind+ind_adjuster;
                
             if my_grid_ind(i) ~= 1
                 if d.dim==2
                     data_p1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2)), pMATLAB.tag, pMATLAB.comm);
                 elseif d.dim==3
                     data_p1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2), send_proc_ind(3)), pMATLAB.tag, pMATLAB.comm);
                 elseif d.dim==4
                     data_p1{i} = MPI_Recv(proc_grid(send_proc_ind(1), send_proc_ind(2), send_proc_ind(3), send_proc_ind(4)), pMATLAB.tag, pMATLAB.comm);
                 end
             end

%now perform the operation - in future will look something like
% if op = '+' (add)
% elseif op = 'mult' (multiply) etc
             if my_grid_ind(i) ~=1    
                 subs = cell(1,d.dim);
                 for s = 1:d.dim
                     subs{s} = ':';
                 end
                 subs{i} = 1:p.overlap(i);
                 d.local(subs{:}) = d.local(subs{:}) + data_p1{i};
             end
  
             if my_grid_ind(i) ~=grid_dims(i)
                 subs = cell(1,d.dim);
                 for s = 1:d.dim
                     subs{s} = ':';
                 end
                 local_dims = size(d.local);
                 subs{i} = (local_dims(i)-p.overlap(i)+1):local_dims(i);
                 d.local(subs{:}) = d.local(subs{:}) + data_m1{i};
             end   
        end %overlap is greater than 0
     end %overlap description defined
 end % inmap

end %synchronize dimension loop

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
