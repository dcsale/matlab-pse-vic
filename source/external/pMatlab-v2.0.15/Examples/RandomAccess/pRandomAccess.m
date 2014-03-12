%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUPS (Giga UPdates per Second) is a measurement that profiles the memory
% architecture of a system and is a measure of performance similar to MFLOPS.
% GUPS is calculated by identifying the number of memory locations that can be
% randomly updated in one second, divided by 1 billion (1e9). The term
% "randomly" means that there is little relationship between one address to be 
% updated and the next, except that they occur in the space of one half the
% total system memory.  An update is a read-modify-write operation on a table
% of 64-bit words.  An address is generated, the value at that address read from
% memory, modified by an integer operation (add, and, or, xor) with a literal
% value, and that new value is written back to memory.
%
% Parameters:
%   * PARALLEL - Turn on distributed arrays.
%
%   * VALIDATE - Enable validation of results.
%
%   * Nup - Number of updates to the table.
%
%   * ErrorRate - Allowed error rate.
%
%   * lgN - Used to set size of table, where N = 2^lgN
%
%   * lgNb - Used to set block size of generating random references
%      where Nb = 2^lgNb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run in serial without distributed arrays, set         
%   PARALLEL = 0
% At the Matlab prompt type
%   pRandomAccess
% To run in serial with distributed arrays, set         
%   PARALLEL = 1
% At the Matlab prompt type
%   pRandomAccess
% To run in parallel with distributed arrays
% at the Matlab prompt type 
%   eval(pRUN('pRandomAccess',2,{}))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment if you want all numbers printed out in hex.
% format hex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS (OK to change)
PARALLEL = 1;  % Turn distributed arrays on or off.
VALIDATE = 1;  % Turn verification on or off.
ErrorRate = 0.01;  % Set allowed error rate.

% Log size of main table (suggested: half of global memory).
%lgN = 27+log2(Np);   % Full memory.
%lgN = 20;   % Performance.
lgN = 25;   % Book
lgN = 18;    % Debug.

% Since the main table size is large, updates are performed in blocks instead
% of one at a time.
% Set log block size.  The official HPC Challenge benchmark sets the block size
% to 1024, i.e. sets the lgNb to 10
lgNb = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETUP

N = 2^lgN;  % Size of update table X.

% Number of updates to table.
% Must run either 4x number of table OR at least 1/4 time to run the Top500
% benchmark (~1000 seconds?)
%Nup = (4 * N);
%Nup = N/4;
Nup = N/16;
%Nup = N/2^(lgN-10-8);

Nb = 2^lgNb;  % Size of update blocks.
Nblocks = Nup / Nb;  % Number of update blocks.

% Create mask that selects low bits we want to use for indexing large table X.
mask = uint64(N-1);

% Create maps.
Xmap = 1;
if PARALLEL
  Xmap = map([1 Np],{},0:Np-1);  % Map for table.
end

tic;
  X = zeros(1,N,Xmap);  % Allocate main table.
  Xloc = uint64(global_ind(X,2) - 1);   % Initialize local part indices.
Talloc = toc;
disp(['Allocation Time (sec)              = ',num2str(Talloc)]);

myX = global_block_range(X,2);  % Get local index range.
allX = global_block_ranges(X,2); % Get all index ranges.

disp(['Distributed table size             = 2^' ...
      num2str(lgN) ' = ' num2str(N) ' words']);
disp(['Distributed table size (bytes)     = ' num2str(N.*8)]);
disp(['Local table size (bytes)           = ' num2str(numel(Xloc).*8)]);
disp(['Number of updates                  = ' num2str(Nup)]);
disp(['Block size (should be 1024)        = ' num2str(Nb)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN BENCHMARK

% Distribute update block indices across processors.
myBLOCK = global_ind(zeros(1, Nblocks,Xmap),2);

% Create a block of starting locations in the random sequence.
ranStarts = RandomAccessStarts( ...
   (myBLOCK(1)-1)*Nb + (0:(Nb-1))*length(myBLOCK) );

% Synchronize start.
tic;
  sync = agg(zeros(1, Np, Xmap));
Tlaunch = toc;
disp(['Launch Time (sec)                  = ',num2str(Tlaunch)]);

% Cache VALIDATE flag.
tempVALIDATE = VALIDATE;
VALIDATE = 0;

% Run core of parallel RandomAccess benchmark, without validation
tic;

%   pRandomAccessSpray;
   pRandomAccessTree;

Trun = toc;
disp(['Run time (sec)                     = ',num2str(Trun)]);

% Compute GUPS.
GUPS = Nup / Trun / 1.e9;
disp(['Giga Updates Per Sec               = ',num2str(GUPS)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATION (in serial or "safe" mode; optional)

% Put VALIDATE flag back.
VALIDATE = tempVALIDATE;

if VALIDATE
  disp('Validating results...');

  tic;

    % Run core of parallel RandomAccess benchmark, with validation
    pRandomAccessSpray;
%    pRandomAccessTree;

    Xloc0 = uint64(global_ind(X,2) - 1);  % Compute errors.
    Nerrors = length(find(Xloc ~= Xloc0));

Tvalidate = toc;
disp(['Validate time (sec)                = ',num2str(Tvalidate)]);

  % Aggregate data back on leader.
  Etable = zeros(1, Np, Xmap);
  Etable = put_local(Etable,Nerrors);
  EtableAll = agg(Etable);

  if (Pid == 0)
    % Total errors.
    erate = sum(EtableAll)/N;
disp(['Error rate                         = ',num2str(erate)]);
    if (erate > ErrorRate)
      disp('Validation Failed')
    else
      disp('Validation Passed')
    end
  end
end
