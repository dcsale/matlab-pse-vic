%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements a simple parallel stream benchmark.  It times the
% following operations:
%   COPY          C(i) = A(i)
%   SCALE         B(i) = q*C(i)
%   ADD           C(i) = A(i) + B(i)
%   TRIAD         A(i) = B(i) + q*C(i)
%
% For parallel systems, we impose the constraint that the size of A, B, and C
% should take up at least half of the total system memory.
%
% Parameters:
%   * PARALLEL - Enable the pMatlab library.
%
%   * N - Length of vectors A, B, and C in number of elements.  Can set this
%      value directly or use lgN.
%
%   * lgN - Used to set N, where N = 2^lgN.
%
%   * NTRIALS - Number of trials to perform for each operation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pStream
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pStream
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pStream',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS (OK to change)

% Turn parallelism on or off
PARALLEL = 1;  % Can be 1 or 0.

% Set number of trials to perform for each operation
NTRIALS = 10;

% Scale data size by number of cpus size
lgN = 25;
lgN = 20; % Debug.
N = 2.^lgN;
%N = floor(1.6 .* 2.^lgN);  % Largest on 1 TX-2500 node.
%N = Np .* 2.^lgN;
%N = Np .* floor(1.6 .* 2.^lgN);  % Largest on Np TX-2500 nodes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP

% Pick initial values.
A0 = 1.0; B0 =2.0; C0 = 0.0;  q = 3.14;

% Compute final values.
ANm1 = ((2.*q + q.^2).^(NTRIALS-1)).*A0;
AN = ((2.*q + q.^2).^(NTRIALS)).*A0;
BN = q.*ANm1;
CN = (1+q).*ANm1;

% Create maps
ABCmap = 1;
SyncMap = 1;
if (PARALLEL)
    % Create map.
    ABCmap = map([1 Np], {},0:Np-1);
    SyncMap = map([Np 1], {},0:Np-1);
end

% Allocate data structures
tic;
    Aloc = local(zeros(1, N, ABCmap)) + A0;
    Bloc = local(zeros(1, N, ABCmap)) + B0;
    Cloc = local(zeros(1, N, ABCmap)) + C0;
Talloc = toc;
disp(['Allocation Time (sec)              = ',num2str(Talloc)]);
Nloc = numel(Cloc);

% Perform barrier synchronization with agg().
%!!!!! DO WE NEED THIS?
tic;
  sync = agg(zeros(1, Np, SyncMap));
Tlaunch = toc;
disp(['Launch Time (sec)                  = ',num2str(Tlaunch)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN PROGRAM
TsumCopy=0.0; TsumScale =0.0; TsumAdd=0.0; TsumTriad=0.0;

for i_trial = 1:NTRIALS

    % COPY
    tic;
        Cloc(:,:) = Aloc;
    TsumCopy = TsumCopy + toc;

    % SCALE
    tic;
        Bloc(:,:) = q*Cloc;
    TsumScale = TsumScale + toc;

    % ADD
    tic;
        Cloc(:,:) = Aloc + Bloc;
    TsumAdd = TsumAdd + toc;

    % TRIAD
    tic;
        Aloc(:,:) = Bloc + q*Cloc;
    TsumTriad = TsumTriad + toc;

end

% Check results.
Aerr = max(abs(Aloc-AN));
Berr = max(abs(Bloc-BN));
Cerr = max(abs(Cloc-CN));


% Compute performance.
copy_BW = 16.*Nloc.*NTRIALS ./ TsumCopy;
scale_BW = 16.*Nloc.*NTRIALS ./ TsumScale;
add_BW = 24.*Nloc.*NTRIALS ./ TsumAdd;
triad_BW = 24.*Nloc.*NTRIALS ./ TsumTriad;

% DISPLAY RESULTS
disp(['A, B, C errors                     = ' num2str([Aerr Berr Cerr])]);
disp(['Np                                 = ' num2str(Np)]);
disp(['Pid                                = ' num2str(Pid)]);
disp(['Global Array size (elem)           = ' num2str(N)]);
disp(['Global Array size (bytes)          = ' num2str(N*8)]);
disp(['Global memory required (bytes)     = ' num2str(N*24)]);
disp(['Local Array size (elem)            = ' num2str(Nloc)]);
disp(['Local Array size (bytes)           = ' num2str(Nloc*8)]);
disp(['Local memory required (bytes)      = ' num2str(Nloc*24)]);
disp(['Number of trials                   = ' num2str(NTRIALS)]);
disp(['Local Copy Bandwidth (MB/sec)      = ' num2str(copy_BW/1e6)]);
disp(['Local Scale Bandwidth (MB/sec)     = ' num2str(scale_BW/1e6)]);
disp(['Local Add Bandwidth (MB/sec)       = ' num2str(add_BW/1e6)]);
disp(['Local Triad Bandwidth (MB/sec)     = ' num2str(triad_BW/1e6)]);
disp(['Global Copy Bandwidth (MB/sec)     ~ ' num2str(Np.*copy_BW/1e6)]);
disp(['Global Scale Bandwidth (MB/sec)    ~ ' num2str(Np.*scale_BW/1e6)]);
disp(['Global Add Bandwidth (MB/sec)      ~ ' num2str(Np.*add_BW/1e6)]);
disp(['Global Triad Bandwidth (MB/sec)    ~ ' num2str(Np.*triad_BW/1e6)]);
