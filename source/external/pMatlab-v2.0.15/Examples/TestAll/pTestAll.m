%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test all the examples.
% AddOne
% Mandelbrot
% ZoomImage
% IO
% Beamformer
% Blurimage
% Stream
% RandomAccess
% FFT
% HPL
% Speedtest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Will need to replace 'clear' with 'clearvars('-except','pMATLAB');'


disp('pAddOne'); cd('../AddOne'); pAddOne; clear('-x','pMATLAB'); pause(2);
disp('pMandelBrot'); cd('../Mandelbrot'); pMandelbrot; clear('-x','pMATLAB'); pause(2);
disp('pZoomImage'); cd('../ZoomImage'); pZoomImage; clear('-x','pMATLAB'); pause(2);
disp('pIO'); cd('../IO'); pIO; clear('-x','pMATLAB'); pause(2);
disp('pBeamformer'); cd('../Beamformer'); pBeamformer; clear('-x','pMATLAB'); pause(2);
disp('pBlurimage'); cd('../Blurimage'); pBlurimage; clear('-x','pMATLAB'); pause(2);
disp('pStream'); cd('../Stream'); pStream; clear('-x','pMATLAB'); pause(2);
disp('pRandomAccess'); cd('../RandomAccess'); pRandomAccess; clear('-x','pMATLAB'); pause(2);
disp('pFFT'); cd('../FFT'); pFFT; clear('-x','pMATLAB'); pause(2);
disp('pHPL'); cd('../HPL'); pHPL; clear('-x','pMATLAB'); pause(2);
if (Np == 2)
  disp('pSpeedtest'); cd('../Speedtest'); pSpeedtest; clear('-x','pMATLAB'); pause(2);
end

cd('../TestAll');
