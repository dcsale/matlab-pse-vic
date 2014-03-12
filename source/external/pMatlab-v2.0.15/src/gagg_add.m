function m = gagg_add(m, n, Ops)
%GAGG_ADD(m) does an operation on arrays based on their types and operator definition
%
%    GAGG(M, N, Destination, Operator)
%      M, N:     Data to be gathered
%      Operator: Operator to be used for the gather operation   
%
% Author: Dr. Chansup Byun
           
%%
name = deblank(class(m));
if strncmp(name, 'Assoc', 5) 
   % Associative class 
   if strncmp(Ops, '+', 1)
      m = m + n;
   elseif strncmp(Ops, '-', 1)
      m = m - n;
   else
      error(['gagg: ' Name ' class does not support the operator, ' Ops]);
   end
elseif strncmp(name, 'double', 6) || strncmp(name, 'single', 6)
   % single or double class
   if strncmp(Ops, 'min', 3)
      m = min(m, n);
   elseif strncmp(Ops, 'max', 3)
      m = max(m, n);
   elseif strncmp(Ops, '+', 1) || strncmp(Ops, 'sum', 3) 
      % Matlab standard plus operation between two same-size arrays
      m = m + n;
   else
      error(['gagg: ' name ' class does not support the operator, ' Ops]);
   end
else
   error(['Current gagg() does not support class ' name ' object.']);
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
