%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pMatlab: Parallel Matlab Toolbox
% Software Engineer: Ms. Nadya Travinin (nt@ll.mit.edu)
% Architect: Dr. Jeremy Kepner (kepner@ll.mit.edu)
% MIT Lincoln Laboratory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN is a generic script for running pMatlab scripts.

% Define number of processors to use.
Ncpus = 4;

% Code to run.
mFile = 'pTestAll';

% Define machines, empty means run locally.
machines = {};

% Define LL Grid machines.

rack_order = [2];

mi = 1;
for j=rack_order
 rack = num2str(j);
%  for i=0:15
  for i=1:15
    is = num2str(i);
    machines(mi) = {['f-' rack '-' is ':/state/partition1/crossmounts/f-' rack '-' is '/kepner']};
    machines_ib(mi) = {['f-' rack '-' is ':/state/partition1/crossmounts/ib-f-' rack '-' is '/kepner']};
    machines_rdma(mi) = {['f-' rack '-' is ':/state/partition1/rdma/f-' rack '-' is '/kepner']};
    mi = mi+1;
  end
end
rack = num2str(rack_order(1));
%machines(1) = {[getenv('HOSTNAME') ':/state/partition1/crossmounts/f-' rack '-0/kepner']};
machines(1) = {[getenv('HOSTNAME') ':/state/partition1/crossmounts/f-' rack '-1/kepner']};
%machines_ib(1) = {[getenv('HOSTNAME') ':/state/partition1/crossmounts/ib-f-' rack '-0/kepner']};
%machines_rdma(1) = {[getenv('HOSTNAME') ':/state/partition1/rdma/f-' rack '-0/kepner']};

%cpus = 'grid';

%cpus = machines_rdma([1 3 4]);
cpus = machines([1 11]);

% Run using local machine and tmp.
%cpus = {[getenv('HOSTNAME') ':/tmp/kepner']};

cpus = {};

% Run the script.
disp(['Running: ' mFile ' on ' num2str(Ncpus) ' cpus']);
eval(pRUN(mFile, Ncpus, cpus));

