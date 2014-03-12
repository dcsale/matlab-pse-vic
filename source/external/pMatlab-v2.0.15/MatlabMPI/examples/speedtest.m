%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script times MPI_Send/MPI_Recv for
% a variety of message sizes.
% To run, start Matlab and type:
%
%   eval( MPI_Run('speedtest',2,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('speedtest',2,{'machine1' 'machine2'}) );
%
% Output will be piped into to
%
%   MatMPI/speedtest.0.out
%   MatMPI/speedtest.0.mat
%   MatMPI/speedtest.1.out
%   MatMPI/speedtest.1.mat
%   ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the number message sizes.
n_message = 20;

% Set the number of trials at each messages size.
n_trial = 4;

% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Modify common directory from default for better performance.
% comm = MatMPI_Comm_dir(comm,'/tmp');
% comm = MatMPI_Comm_dir(comm,'/wulf/share/kepner');
% comm = MatMPI_Comm_dir(comm,'/gigabit/node-a');

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Do a synchronized start.
% starter_rank = 0;
% delay = 30;  % Seconds
% synch_start(comm,starter_rank,delay);

if(comm_size < 2)
 disp('ERROR: too few processors (need at least 2)');
 exit;
end

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Set who is source and who is destination.
source = my_rank - 1;
if (source < 0)
  source = comm_size - 1;
end
dest = my_rank + 1;
if (dest >= comm_size)
  dest = 0;
end


% Create a unique tag id for this message (very important in Matlab MPI!).
tag = 1;

% Create timing matrices.
start_time = zeros(n_trial,n_message);
end_time = start_time;

% Compute message sizes.
p = 1:n_message;
message_size = 2.^p;
byte_size = 8.*message_size;


% Get a zero clock.
zero_clock = clock;


% Loop over each message size.
for i_message = 1:n_message


  % Create message.
  send_data = zeros(1,message_size(i_message)) + my_rank;

  for i_trial = 1:n_trial

    % Get start time for this message.
    start_time(i_trial,i_message) = etime(clock,zero_clock);

    % Send data.
    MPI_Send( dest, tag, comm, send_data );

    % Recieve data.
    recv_data = MPI_Recv( source, tag, comm );

    % Get end time for the this message.
    end_time(i_trial,i_message) = etime(clock,zero_clock);

    total_time = end_time(i_trial,i_message) - start_time(i_trial,i_message);

    % Check data.
    if(any(recv_data ~= source))
      disp('ERROR: incorrect data sent.');
      exit;
    end

    % Increment message tag.
    tag = tag + 1;

  end
end

% Compute bandwidth.
total_time = end_time - start_time;
byte_size_matrix = repmat(byte_size,n_trial,1);
bandwidth = 2.*byte_size_matrix./total_time;

% Write data to a file.
outfile = ['speedtest.',num2str(my_rank),'.mat'];
save(outfile,'byte_size','start_time','end_time','total_time','byte_size_matrix','bandwidth');

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

