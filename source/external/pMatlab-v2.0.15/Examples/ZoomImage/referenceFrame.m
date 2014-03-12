function refFrame = referenceFrame(n,a,b);

  % Create axis.
  x_axis = linspace(-1,1,n);
  y_axis = linspace(-1,1,n);
  nx = length(x_axis);
  
  
  % Create x and y.
  [x y] = meshgrid(x_axis,y_axis);

  rt2 = 2.^0.5; 

  % Create Reference frame.
  refFrame = (y > x - rt2*a) & (y < x + rt2*a) ...
          &  (y < rt2*b - x) & (y > -rt2*b - x);

  % Add some noise (shouldn't add noise).
  refFrame = refFrame + 0.1*(rand(size(refFrame)) - 0.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2005, Massachusetts Institute of Technology All rights     %
%reserved.                                                                %
%                                                                         %
%Redistribution and use in source and binary forms, with or without       %
%modification, are permitted provided that the following conditions are   %
%met:                                                                     %
%     * Redistributions of source code must retain the above copyright    %
%       notice, this list of conditions and the following disclaimer.     %
%     * Redistributions in binary form must reproduce the above copyright %
%       notice, this list of conditions and the following disclaimer in   %
%       the documentation and/or other materials provided with the        %
%       distribution.                                                     %
%     * Neither the name of the Massachusetts Institute of Technology nor %
%       the names of its contributors may be used to endorse or promote   %
%       products derived from this software without specific prior written% 
%       permission.                                                       %
%                                                                         %
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS  %
%IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,%
%THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR   %
%PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR         %
%CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,    %
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,      %
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR       %
%PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF   %
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     %
%NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS       %
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

