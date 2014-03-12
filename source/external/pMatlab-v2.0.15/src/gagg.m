function m = gagg(m, varargin)
%GAGG(m) gathers the parts of an array on the leader processor
%    or does global reduction on a variable on the leader processor.
%    GAGG(M) gathers the parts of an array on the remote processes
%    into a whole array on the leader process. 
%    If the current processor is the LEADER, returns 
%    the aggreagated array, otherwise, returns the 
%    local part.
%    This functions increments GLOBAL message TAG.
%
%    GAGG(M)
%    GAGG(M, Destination, Operator)
%      M:  Data to be gathered
%      Destination:  Destination where the data is to be gathered
%                    With only M, Destination is Pid=0 by default
%      Operator:     Operator to be used for the gather operation   
%
% Author: Dr. Chansup Byun
%
% globals
global pMATLAB;

% disp('Debug: Entered GAGG()');   
if isfield(pMATLAB, 'tag_num')
  if Np == 1
     return
  else 
     %increment tag
     pMATLAB.tag_num = pMATLAB.tag_num+1;
     pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));
  end
else
  % make sure it works for serial mode
  % disp('Debug: Exited GAGG() with No Op. Not a pMatlab run.');   
  return;
end

Ops  = '+'; % Default operator
plist = 0:Np-1; % List of process Pids

if nargin == 1 % total number of arguments
   Dest = 0;   % Destination
elseif nargin == 2
   Dest = varargin{1};
elseif nargin == 3
   Dest = varargin{1};
   Ops  = varargin{2};
   %    disp(deblank(class(Dest)));
   %    disp(deblank(class(Ops)));
elseif nargin > 4
   error('gagg: Incorrect number of inputs');
end
name = class(m);
if strncmp(name, 'Assoc', 5) || strncmp(name, 'double', 6) || strncmp(name, 'single', 6)
   % supported
   % disp(['Class name: ' name]);
else
   error(['Current gagg() does not support class ' name ' object.']);
end

%% Gathering based on binary tree shown below
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
pidList = 0:POTN-1; % Store Pids that participates with msg communication at the current level
pidPost = 1:POTN;   % Pid position in pidList
%
btMax = length(tree);  % the max level of binary tree for a given POTN
bt = 1;  % the starting binary tree level

% Walk up the binary tree.
while (bt <= btMax)
   % Find my Pid position in pidList
   % (Search is limited to the active Pid list)
   % (There are ficticious Pid numbers >= Np when Np != POTN)
   % pidList(1:length(pidPost))
   myPidPos = find( pidList(1:length(pidPost)) == Pid);
   if (myPidPos)
      % Pid participates in communication at this level
      if ( mod(myPidPos,2) )
         % Receive message from my right neighbor, pidList(myPidPos+1)
         % only if my right neighbor is real Pid
         fromRank = pidList(myPidPos+1);
         % disp(['  myPidPos+1 = ' num2str(myPidPos+1) ', fromRank = ' num2str(fromRank)]);
         if fromRank < Np % Receive data only if fromRank is a REAL Pid
            % disp(['gagg() recv: Pid = ' num2str(Pid) ' fromRank ' num2str(fromRank)]);
            m = gagg_add(m, RecvMsg(fromRank,pMATLAB.tag), Ops);
         end
      else
         % Send message to my left neighbor, pidList(myPidPos-1)
         toRank = pidList(myPidPos-1);
         % res = whos('m');
         % disp(['gagg() sent: Pid = ' num2str(Pid) ' toRank ' num2str(toRank)]);
         % disp(['             MsgSize(bytes) = ' num2str(res.bytes)]);
         if toRank < Np % Send data only if toRank is a REAL Pid
            SendMsg(toRank, pMATLAB.tag, m)
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
% disp('Debug: Exited GAGG()');
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
