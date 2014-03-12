function mat = agg(d)
% Hierarchical AGG aggregates the parts of a distributed matrix on 
%    the leader processor using an extended binary tree.
%    HAGG(D) is a generalization of BAGG(D) when Np is NOT a power of two
%    This functions increments GLOBAL message TAG.
%
%    HAGG(D) is renamed to be AGG(D)
%   
%    NOTE: Currently, it doesn't matter whether or not the leader is in the
%    map - the global matrix is returned on the leader regardless. 
%
% Author: Dr. Chansup Byun (cbyun@ll.mit.edu)
% Date:   May 3, 2010
          
%globals
global pMATLAB;

%increment tag
pMATLAB.tag_num = pMATLAB.tag_num+1;
pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));

%% Aggregation based on binary tree shown below 
%
%  (Rank 0 will collect the result)
%                          0
%             0                          8               k = 4 (8 Units)
%       0           4            8              12       k = 3 (4 Units)
%    0     2     4     6     8      10      12     14    k = 2 (2 Units)
%  0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  k = 1 (1 Unit)
%
tree = factor(Np);
% Check if Np is power of two
if (size(find(tree==2))==size(tree))
   POTN = Np;
   % disp(['The Np is ', num2str(Np)]);
else
   % Np is NOT a power of two
   % Find the next power of two number
   k = 0;
   POTN = 1;
   while ( Np > POTN)
      k = k + 1;
      POTN = POTN * 2;
   end
   % disp(['The next power of two number is ', num2str(POTN)]);
   % Generate a new binary tree using the next power of two number
   tree = factor(POTN);
end
pidList = 0:POTN-1; % Store Pid's that participates with msg communication at the current level
pidPost = 1:POTN;   % Pid position in pidList
% 
btMax = length(tree);  % the max level of binary tree for a given POTN
bt = 1;  % the starting binary tree level

if d.dim > 4
   error('DMAT/AGG: Only up to 4-D objects currently supported');
end

dim = size(d.map.grid);
totalProcs = 1;
for i = 1:length(dim)
    totalProcs = totalProcs * dim(i);
end
if (totalProcs < Np)
    % If totalProcs doesn't match with Np, use the old serialized agg()
    mat = oagg(d);
    return
end
%dim(1) - number of grid rows, dim(2) - number of grid col
temp_mat = cell(dim); % If commented, gets the following error w/ Octave.
%                       error: invalid assignment to cs-list outside multiple assignment.
%                       error: assignment to cell array failed
%                       error: assignment failed, or no method for `cell = matrix'

% Generate relation between process rank and grid map
gridIndex = mapGridRank(d);
   
if (pMATLAB.my_rank == pMATLAB.leader) % hagg() leader
   % local data
   if d.dim==2
      % Two dimensional array
      %
      [i j] = findGridIndex(d.dim,pMATLAB.leader,gridIndex);
      temp_mat{i,j} = d.local;
   elseif d.dim==3
      % Three dimensional array
      %
      [i j k] = findGridIndex(d.dim,pMATLAB.leader,gridIndex);
      % disp(['hagg() leader: i j k = ' num2str(i) ', ' num2str(j) ...
      %       ', ' num2str(k)]);
      temp_mat{i,j,k} = d.local;
   elseif d.dim==4
      % Four dimensional array
      %
      [i j k m] = findGridIndex(d.dim,pMATLAB.leader,gridIndex);
      temp_mat{i,j,k,m} = d.local;
   end
else
   % Non-leader returns the local data
   mat = d.local;
   % Non-leader prepares to send its local data
   sendBuf{1,1} = d.local;
   imsgLast = 0;  % pointer to manage send buffer location
end
   
if (Np == 1) 
    mat = d.local;
    return;
end

% Walk up the binary tree.
while (bt <= btMax)
   % Compute msg units transferred at this level
   msgUnit = 2^(bt-1);
   % Find my Pid position in pidList
   % (Search is limited to the active Pid list)
   % (There are ficticious Pid numbers >= Np when Np != POTN)
   % pidList(1:length(pidPost))
   myPidPos = find( pidList(1:length(pidPost)) == Pid);
   if (myPidPos)
      % Pid participates in communication at this level
      if ( mod(myPidPos,2) )
         % Receive message from my right neighbor, pidList(myPidPos+1)
         fromRank = pidList(myPidPos+1);
         % disp(['  myPidPos+1 = ' num2str(myPidPos+1) ', fromRank = ' num2str(fromRank)]);
         if inmap(d.map, fromRank) % Only receive data if fromRank is in the map
            % disp(['hagg() recv: Pid = ' num2str(Pid) ', fromRank ' ...
            %      num2str(fromRank) ' w/ msg unit = ' num2str(msgUnit)]);
            recvBuf = MPI_Recv(fromRank, pMATLAB.tag, pMATLAB.comm);
            % disp(['  size(recvBuf,2) = ' num2str(size(recvBuf,2))]);
            if (pMATLAB.my_rank == pMATLAB.leader) 
               % hagg() leader puts the received msg into temp dmat
               for imsg=1:size(recvBuf,2)
                   if d.dim==2
                      % Two dimensional array
                      % disp(['fromRank+imsg-1 = ' num2str(fromRank+imsg-1)]);
                      [i j] = findGridIndex(d.dim,fromRank+imsg-1,gridIndex);
                      temp_mat{i,j} = recvBuf{1,imsg};
                   elseif d.dim==3
                      % Three dimensional array
                      [i j k] = findGridIndex(d.dim,fromRank+imsg-1,gridIndex);
                      temp_mat{i,j,k} = recvBuf{1,imsg};
                   elseif d.dim==4
                      % Four dimensional array
                      [i j k m] = findGridIndex(d.dim,fromRank+imsg-1,gridIndex);
                      % disp(['i j k m = ' num2str(i), ', ', num2str(j) ...
                      %       ', ' num2str(k) ', ' num2str(m)]);
                      temp_mat{i,j,k,m} = recvBuf{1,imsg};
                   end
               end
            else
               % Others puts the received msg to a buffer for the next level
               % msg aggregation. Always store the local first and then, the
               % msg from the right neighbors sequentially
               for imsg=1:size(recvBuf,2)
                  sendBuf{1,imsgLast + imsg+1} = recvBuf{1,imsg};
               end
               imsgLast = imsgLast + size(recvBuf,2);
            end 
         end
      else
         % Send message to my left neighbor, pidList(myPidPos-1)
         toRank = pidList(myPidPos-1);
         if inmap(d.map, pMATLAB.my_rank) % Only send data if processor is in the map
            % disp(['hagg() sent: Pid = ' num2str(Pid) ', toRank ' ...
            %       num2str(toRank) ' w/ msg unit = ' num2str(msgUnit)]);
            % disp(['  size(sendBuf,2) = ' num2str(size(sendBuf,2))]);
            MPI_Send(toRank, pMATLAB.tag, pMATLAB.comm, sendBuf);
         end
      end
   end
   %
   % Prepare for the next level (reduce size in half)
   bt = bt + 1;
   pidPost = 1:length(pidPost)/2;
   for in = 1:length(pidPost)
       pidList(in) = pidList(2*in-1);
   end
end
% Finalize on the leader processor
if (pMATLAB.my_rank == pMATLAB.leader) % hagg() leader
   % Reconstruct the matrix from the local pieces
   % This is a NO-OP for block distributions 
   mat = reconstruct(d.pitfalls,  d.map.grid, temp_mat, d.size); 
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
