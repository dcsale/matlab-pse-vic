% Script that calls all the unit tests.
% To run:
%
%  (1) Make sure your settings in MatMPI_Comm_settings.m
%      are correct for your system.
%
%  (2) Make sure the machines listed in machine.m are
%      are correct for your system
%
%  (3) Launch Matlab and type 'unit_test_all'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERY VERY VERY IMPORTANT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the HOST is the first machine in your
% machines list, you must run it in the
% background, i.e.
%   machines = {'host&' 'machine2' ... }

% Change list of cpus to suit your environment.
% cpus = {'SLAVE&:pc'};
% cpus = {'SLAVE&:pc:C:\Documents and Settings\kepner\tmp'}

% SunOS
% cpus = {'hagar&','buddy'};

% Change list of cpus to suit your environment.
%cpus = {'node-a&:/gigabit/node-a/kepner' ...
%        'node-b:/gigabit/node-b/kepner' ...
%        'node-c:/gigabit/node-c/kepner' ...
%        'node-d:/gigabit/node-d/kepner' ...
%        'node-e:/gigabit/node-e/kepner' ...
%        'node-f:/gigabit/node-f/kepner' ...
%        'node-g:/gigabit/node-g/kepner' ...
%        'node-h:/gigabit/node-h/kepner'};


% Abort left over jobs.
MPI_Abort;
pause(2.0);

% Delete left over MPI directory
MatMPI_Delete_all;
pause(2.0);

% Set initial status.
all_status = 'SUCCESS';

% Unit test all the scripts.
status = unit_test( 'broadcast', 4, cpus, 200, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'basic_app', 4, cpus, 400, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'basic_app2', 4, cpus, 400, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'basic_app3', 4, cpus, 400, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

%status = unit_test( 'basic_app4', 4, cpus, 400, 20 );
%if (strcmp(status,'FAIL'))
%  all_status = 'FAIL';
%end

status = unit_test( 'xbasic', 4, cpus, 200, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'basic', 2, cpus, 200, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'multi_basic', 2, cpus, 200, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'probe', 2, cpus, 200, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end


status = unit_test( 'speedtest', 2, cpus, 400, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

status = unit_test( 'blurimage', 4, cpus, 400, 20 );
if (strcmp(status,'FAIL'))
  all_status = 'FAIL';
end

disp(['All tests: ' all_status]);

