function m = gagg_all(m, varargin)
%GAGG_ALL(m) broadcasts the gagg() results to all the processors.
%    GAGG_ALL(M) broadcasts the gagg() results to all the processors.
%
%    This functions increments GLOBAL message TAG.
%
%    GAGG_ALL(M)
%    GAGG_ALL(M, SourceOfBroadcast)
%    GAGG_ALL(M, SourceOfBroadcast, Operator)
%      M:  Data to be gathered & broadcasted
%      SourceOfBroadcast:  Source of broadcasting where the data 
%                    are to be gathered and broadcasted.
%                    With only M, SourceOfBroadcast is Pid=0 by default
%      Operator:     Operator to be used for the gather operation   
%
% Author: Dr. Chansup Byun

%globals
global pMATLAB;
if isfield(pMATLAB, 'tag_num')
  %increment tag
  pMATLAB.tag_num = pMATLAB.tag_num+1;
  pMATLAB.tag = strcat('tag-', num2str(pMATLAB.tag_num));
else
  return
end

Ops  = '+'; % Default operator
if nargin == 1 % total number of arguments
   srcBC = 0;   % SourceOfBroadcast
elseif nargin == 2
   srcBC = varargin{1};
elseif nargin == 3
   srcBC = varargin{1};
   Ops  = varargin{2};
   %    disp(deblank(class(srcBC)));
   %    disp(deblank(class(Ops)));
elseif nargin > 4
   error('gagg: Incorrect number of inputs');
end
m = gagg(m,srcBC,Ops);
if Np == 1
   return
end
           
%
% Scatter operation
%
% disp('Inside gagg_all(): calling BcastMsg()...');
% whos m;
m = BcastMsg(0,pMATLAB.tag,m);

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
