%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates a distributed matrix. Writes it
% out in a scalable way and the reads it back in.
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pIO
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pIO
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pIO',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the problem dimensions and filename.
N = 2^16;   M = 128;   FILE = './dat/pIOvector';
N = 2^14;  % Debug.
PARALLEL = 1;  % Set control flag.
Xmap = 1;  Ymap = 1;  % Create serial maps.
if (PARALLEL)
  Xmap = map([1 Np],{},0:Np-1);  Ymap = Xmap;  % Create parallel maps.
end
Xrand = rand(N,M,Xmap);    % Create distributed array.
Yrand = zeros(N,M,Ymap);  % Create distributed array.
tic;   % Start clock.
WriteMatrixParallel(Xrand,FILE);   % Save files.
Twrite = toc;   % Stop clock.
disp(['Write Time (sec)                   = ',num2str(Twrite)]);
tic;   % Start clock.
Yrand = ReadMatrixParallel(Yrand,FILE);  % Read files.
Tread = toc;    % Stop clock.
disp(['Read Time (sec)                    = ',num2str(Tread)]);

% Compare results.
max_difference = max(max(abs( local(Xrand) - local(Yrand) )));

if (max_difference > 0)
  disp('ERROR');
  max_difference
else
%  disp('SUCCESS');
end
