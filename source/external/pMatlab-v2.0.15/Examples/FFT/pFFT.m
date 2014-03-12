%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements a simple 1D FFT benchmark using a 2D approach, which is
% the most common way of implementing a 1D FFT in parallel.  This benchmark
% contains a parallel implementation of this algorithm.
%
% Parameters:
%   * PARALLEL - Enable the pMatlab library.
%
%   * VALIDATE - Enable validation of FFT result.
%
%   * ERROR_LIMIT - Error limit used in validation.
%
%   * N - length of vector to FFT, must be divisible by Np^2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pFFT
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pFFT
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pFFT',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARALLEL=1;                % Turn parallelism on and off.
N = 2^25*Np;  M = N/Np;       % Set vector/matrix dimensions.
N = 2^20*Np;  M = N/Np;    % Debug.
VALIDATE = 0;              % Turn validation on or off.
ErrorRate = eps;           % Set error to machine precision.

Xmap = 1;                  % Serial map.
if PARALLEL
  Xmap = map([1 Np],{},0:Np-1);    % Parallel map.
end

disp(['Np                                 = ' num2str(Np)]);
disp(['Pid                                = ' num2str(Pid)]);
disp(['Distributed vector size (words)    = ' num2str(N)]);
disp(['Distributed vector size (bytes)    = ' num2str(N.*16)]);

tic;
  X = complex(rand(1,N,Xmap),rand(1,N,Xmap));    % Create distributed vector.
  Xloc = local(X);                     % Get local part.
  Xshell = put_local(X,0);             % Create an empty X.

  if VALIDATE
    Xloc(:,:) = complex(0,0);                % Reset to zero.
    phases = floor([3 log2(N) sqrt(N)]);     % Create wave phases.
    for i = 1:length(phases)
      omega = (2*pi.*phases(i)./N).* global_ind(X,2);   % Compute angle.
      Xloc = Xloc + complex(cos(omega),sin(omega));     % Add wave to Xloc.
    end
    X = put_local(X,Xloc);               % Insert back into X.
  end

  % Create twiddle factors
  omega = (-2*pi*Pid/N) .* (0:M-1);
  omega = complex(cos(omega),sin(omega));
Talloc = toc;

disp(['Local vector size (bytes)          = ' num2str(numel(Xloc).*16)]);
disp(['Allocation Time (sec)              = ',num2str(Talloc)]);

tic;
  sync = agg(zeros(1, Np, Xmap));      % Synchronize start.
Tlaunch = toc;
disp(['Launch Time (sec)                  = ',num2str(Tlaunch)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN BENCHMARK

tic;
  Xloc = reshape(Xloc,Np,M/Np);            % Reshape local part into a matrix.
  X = put_local(zeros(Np,M,Xmap),Xloc);
Tcomp = toc;
disp('Begin 1st CornerTurn');
tic;
  X = transpose_grid(X);        % Redistribute along 1st dimension.
  Xloc = local(X);
Tcomm = toc;
disp('Begin FFT of 2nd Dimension');
tic;
  Xloc = fft(Xloc,[],2);        % FFT 2nd dimension.
  Xloc = omega .* Xloc;         % Multiply by twiddle factors.
  X = put_local(X,Xloc);
Tcomp = Tcomp + toc;
disp('Begin 2nd CornerTurn');
tic;
  X = transpose_grid(X);        % Redistribute along 2nd dimension.
  Xloc = local(X);
Tcomm = Tcomm + toc;
tic;
  Xloc = fft(Xloc,[],1);        % FFT 1st dimension.
  X = put_local(X,Xloc);
Tcomp = Tcomp + toc;
tic;
  X = transpose_grid(X);        % Redistribute along 1st dimension.
  Xloc = local(X);
  X = put_local(Xshell,Xloc);   % Insert back into vector.
Tcomm = Tcomm + toc;
Trun = Tcomp + Tcomm; 

GigaFlops = 5*N*log2(N)/Trun/1.e9;   % Performance in gigaflops 
GigaByteSec = 3*16*N/Tcomm/1.e9;  % Communication bandwidth.

disp(['Compute time (sec)                 = ',num2str(Tcomp)]);
disp(['Communication time (sec)           = ',num2str(Tcomm)]);
disp(['Run time (sec)                     = ',num2str(Trun)]);
disp(['Performance (Gigaflops)            = ',num2str(GigaFlops)]);
disp(['Bandwidth (Gigabytes/sec)          = ',num2str(GigaByteSec)]);

if VALIDATE
  myX = global_block_range(X,2);     % Get global index.
  iphase = phases + 1;               % Get phase index.
  myiphase = iphase( (myX(1) <= iphase) & (iphase <= myX(2)) ) - myX(1) + 1;  
  Errloc = max(abs(Xloc))./(N.^1.5);
  if not(isempty(myiphase))
    Zloc = sparse(myiphase,1,N,M,1);
    Errloc = max(abs(abs(Xloc).' - Zloc))./(N.^1.5);
  end
disp(['Max local error                    = ',num2str(Errloc)]);
  Err = put_local(zeros(1,Np,Xmap),Errloc);
  Erragg = agg(Err);
  if (Pid == 0)
    maxErr = max(Erragg);
disp(['Max error                          = ',num2str(maxErr)]);
    if (maxErr > ErrorRate)
      disp('Validation Failed')
    else
      disp('Validation Passed')
    end
  end
end
