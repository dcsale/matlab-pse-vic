%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Matlab MPI script that
% sends a matrix to another processor.
%
% To run, start Matlab and type:
%
%   eval( MPI_Run('multi_basic',2,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('multi_basic',2,{'machine1' 'machine2'}) );
%
% Output will be piped into two files:
%
%   MatMPI/multi_basic.0.out
%   MatMPI/multi_basic.1.out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set number of tries.
n_tries = 3;

% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Modify common directory from default for better performance.
% comm = MatMPI_Comm_dir(comm,'/tmp');
% comm = MatMPI_Comm_dir(comm,'/gigabit/node-a');

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Set who is source and who is destination.
source = 1;
dest = 0;

% Create a unique tag id for this message (very important in Matlab MPI!).
tag = 1;

% Get a zero clock.
zero_clock = clock;

% Check to make sure that we have exactly two processors.
if(comm_size == 2)
  for i=1:n_tries
    % Compute message tag.
    i_tag = tag+i;

    % Get start time for this message.
    start_time = etime(clock,zero_clock);

    % Source.
    if (my_rank == source)
      % Create data.
      data = 1:10;
      % Send it.
      MPI_Send( dest, i_tag, comm, data, data );
    end
    % Destination.
    if (my_rank == dest)
      % Receive data.
      [data data1] = MPI_Recv( source, i_tag, comm );
      % Check data.
      if(any((data  - (1:10)) ~= 0))
        disp('ERROR: incorrect data sent.');
        exit;
      end
    end

    % Get end time for the this message.
    end_time = etime(clock,zero_clock);

    % Print total time.
    total_time = end_time - start_time

  end
end

% Finalize Matlab MPI.
disp('SUCCESS');
MPI_Finalize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2002 Massachusetts Institute of Technology
% 
% Permission is herby granted, without payment, to copy, modify, display
% and distribute this software and its documentation, if any, for any
% purpose, provided that the above copyright notices and the following
% three paragraphs appear in all copies of this software.  Use of this
% software constitutes acceptance of these terms and conditions.
%
% IN NO EVENT SHALL MIT BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF MIT HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
% MIT SPECIFICALLY DISCLAIMS ANY EXPRESS OR IMPLIED WARRANTIES INCLUDING,
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
% FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.
%
% THIS SOFTWARE IS PROVIDED "AS IS," MIT HAS NO OBLIGATION TO PROVIDE
% MAINTENANCE, SUPPORT, UPDATE, ENHANCEMENTS, OR MODIFICATIONS.

