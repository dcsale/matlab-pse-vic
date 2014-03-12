%pMATLAB_Init.m
%Initializes all the needed MPI variables, such as
%   number of processors, current processor's rank and
%   which processor is the leader. All of the variables 
%   necessary for communication are stored in the pMATLAB 
%   structure.
%  
%   Fields of the pMATLAB structure:
%       comm        - contains the MatlabMPI communicator
%       comm_size   - size of communicator, i.e. number of processors
%       my_rank     - rank of the local processor
%       leader      - processor with rank 0
%       pList       - list of participating processors
%       tag         - message tag
%       tag_num     - number of messages sent; synchronized accross all
%                     processors in pList
%   Fields to be potentially added in the future  
%       num_tasks   - number of tasks (scopes) created from the beginning of the program 
%       curr_task   - current task (scope)
%       scopes      - contains a cell array of communication scopes; each entry in
%       the array is a struct with the current fields of the pMATLAB
%       structure plus the task_num field
%
% Author:   Nadya Travinin

%global variables
global pMATLAB MPI_COMM_WORLD;

%Initialize MatlabMPI
MPI_Init;

%initialize pMATLAB struct
pMATLAB.comm = MPI_COMM_WORLD; %MPI communicator

pMATLAB.comm_size = MPI_Comm_size(pMATLAB.comm); %comm size
pMATLAB.my_rank = MPI_Comm_rank(pMATLAB.comm); %my processor rank
pMATLAB.leader=0; %set the leader
pMATLAB.host = MatMPI_Host_rank(pMATLAB.comm); %stores the rank of the interactive process; if all processes are running in the background, returns -1 
pMATLAB.pList = [0:pMATLAB.comm_size-1]; %set processor list
pMATLAB.tag_num = 0;
%message tag - MUST be unique for each message
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

%Uncomment if you want to save the messages that were sent.
%pMATLAB.comm = MatMPI_Save_messages(pMATLAB.comm,1);

% %***DATA STRUCTURE NEEDED FOR TASK SCOPING - currently not used***
% pMATLAB.scopes = {};
% pMATLAB.num_tasks = 1; %initially, just the global task (scope)
% pMATLAB.curr_task = 1; %current task (scope)
% global_scope.comm_size = MPI_Comm_size(pMATLAB.comm);
% global_scope.my_rank = MPI_Comm_rank(pMATLAB.comm); 
% global_scope.leader = 0;
% global_scope.pList = [0:global_scope.comm_size-1];
% global_scope.task_num = pMATLAB.num_tasks;
% global_scope.tag_num = 0;
% global_scope.tag = strcat('scope', num2str(global_scope.task_num), '-', num2str(global_scope.tag_num));
% pMATLAB.scopes(num_tasks) = {global_scope};
% %*****************************************************************

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