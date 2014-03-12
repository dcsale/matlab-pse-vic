%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of the HPC Challenge Higher Performance Linpack
% benchmark which solve the equation Ax = b.
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pHPL
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pHPL
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pHPL',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn parallelism on or off.
PARALLEL = 1;  % Can be 1 or 0.
VERIFY = 1;

Amap = 1;                          % Serial map.
if PARALLEL
  Amap = map([1 Np],{},0:Np-1);    % Parallel map.
end

%N = 2^14;  % Largest on a single beta grid node ~600 seconds on 1 cpu.
%N = Np*floor(floor(10000*sqrt(Np))/Np)
%N = Np*floor(floor(12500*sqrt(Np))/Np)
%N = 1500;
%N = 2^12;
%N = 2^10;
%N = Np*floor(floor((2^12)*sqrt(Np))/Np);
N = Np*floor(floor(5000*sqrt(Np))/Np);
N = 240;  % Debug.


disp(['Np                                 = ' num2str(Np)]);
disp(['Pid                                = ' num2str(Pid)]);
disp(['Distributed matrix size            = ' num2str(N) '^2 words']);
disp(['Distributed matrix size (bytes)    = ' num2str(N.^2.*8)]);

tic;
  b = rand(N,1) - 0.5;        % Create a replicated column vector.
  A = rand(N,N,Amap) - 0.5;   % Create a distributed matrix.
Talloc = toc;

disp(['Local matrix size (bytes)          = ' num2str(numel(local(A)).*8)]);
disp(['Allocation Time (sec)              = ',num2str(Talloc)]);

tic;
  sync = agg(zeros(1, Np, Amap));      % Synchronize start.
Tlaunch = toc;
disp(['Launch Time (sec)                  = ',num2str(Tlaunch)]);

tic
  x = pLUsolve(A,b);                   % Solve  A x = b.
Trun = toc;
GigaFlops = (2/3*N^3 + 3/2*N^2)/Trun/1.e9;  % Performance in gigaflops.

disp(['Run time (sec)                     = ',num2str(Trun)]);
disp(['Performance (Gigaflops)            = ',num2str(GigaFlops)]);

if (VERIFY)
  % Scaled residuals
  A = agg(A);
  if (Pid == 0)
    r=A*x - b;

    normA1 = max(sum(abs(local(A))));
    normAInf = max(sum(abs(local(A)),2));

    %res0 = norm(r, Inf);
    res0 = max(sum(abs(local(r)),2));
    %res1 = res0 / (eps * norm(A, 1) * N);
    res1 = res0 / (eps * normA1 * N);
    %res2 = res0 / (eps * norm(A, 1) * norm(x, 1));
    res2 = res0 / (eps * normA1 * norm(x, 1));
    %res3 = res0 / (eps * norm(A, Inf) * norm(x, Inf) * N);
    res3 = res0 / (eps * normAInf * norm(x, Inf) * N);

    threshold = 16;

    if max([res1, res2, res3]) < threshold
     disp('Verification Passed');
    else
     disp('Failure');
    end
  end
end
