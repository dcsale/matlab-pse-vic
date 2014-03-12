function zoomedFrames = zoomFrames(refFrame,scaleFactor,blurSigma)
% ZOOMFRAMES:  Zooms in on a reference frame.
%    Usage: zoomedFrames =zoomFrames(refFrame,scaleFactor,blurSigma)
%
%    Input parameters are as follows:
%       refFrames  = reference image.
%       scaleFactor= vector of zoome scale factors.
%       blurSigma  = standard deviation of blur kernel (in pixels)
%
%    Output variables are as follows:
%       zoomeFrames = 3D array of output images.
%

% Linear dimension of cropping window  (pixels)
windowSize = size(refFrame);
windowSize = windowSize(1);

% Allocate output frames.
zoomedFrames = zeros(windowSize,windowSize,length(scaleFactor));

% Estimate frames at selected ranges.
for i=1:length(scaleFactor)
      
   % Generate point spread function.
   % Pre-selected to be 5-sigma.
   nelem = ceil(scaleFactor(i)*(5*blurSigma));
   if mod(nelem,2) == 0
      % Make odd.
      nelem = nelem+1;
   end
   % Create gaussian kernel.
   x = (0:(nelem-1)) - (nelem-1)./2;
   o = ones(1,nelem);
   r2 = ((x' * o).^2 + (o' * x).^2);
   h = (1./sqrt(2.*pi.*blurSigma)).*exp(r2./blurSigma.^2);
   
   % Calculate blurred image (at requested range).
   blurFrame = conv2(h,refFrame,'full');
      
   % Simulate sampled image at requested range.
   zoomedFrames(:,:,i) = imdown(blurFrame,windowSize,scaleFactor(i));
   
end % End of frame simulation.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function downimg = imdown(I,res,alpha)

% IMDOWN Resamples a grayscale image.
%    DOWNIMG = IMDOWN(I,R,A) Resamples the inpu
%    image I using bilinear interpolation using a sampling
%    lattice with R point samples at A times the original
%    sample spacing.  Extrapolated values are set to zero.

% Create the down-sampling lattice.
npix = size(I,1);
x = alpha*([1:res]-((res+1)/2))+((res+1)/2)+(npix-res)/2;
y = alpha*([1:res]-((res+1)/2))+((res+1)/2)+(npix-res)/2;
[XI,YI] = meshgrid(x,y);

% Down-sample the input grayscale image.
%downimg = interp2(im2double(I),XI,YI,'bilinear');
%downimg = interp2(double(I),XI,YI,'bilinear');
downimg = interp2(double(I),XI,YI,'linear');
downimg(find(isnan(downimg))) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

