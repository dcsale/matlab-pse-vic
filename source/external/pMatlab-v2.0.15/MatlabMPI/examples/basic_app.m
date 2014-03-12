%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic example of a typical application MatlabMPI
% might be used for.  Takes a matrix, breaks it up
% among processors, does a computation, and gathers
% the results.
%
% To run, start Matlab and type:
%
%   eval( MPI_Run('basic_app',4,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('basic_app',4,{'machine1' 'machine2' 'machine3' 'machine4'}) );
%
% Output will be piped into 4 files:
%
%   MatMPI/basic_app.0.out
%   MatMPI/basic_app.1.out
%   MatMPI/basic_app.2.out
%   MatMPI/basic_app.3.out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatlabMPI
% Dr. Jeremy Kepner
% MIT Lincoln Laboratory
% kepner@ll.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MPI.
MPI_Init;

% Create communicator.
comm = MPI_COMM_WORLD;

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Wait momentarily.
pause(2.0);

% Set who is the leader
leader = 0;

% Create base message tags.
coefs_tag = 10000;
input_tag = 20000;
output_tag = 30000;

% Set data sizes.
N1 = 1000;
N2 = 70;

% Leader.
if (my_rank == leader)
  % Create coefficient data.
  coefs = ones(N1,1);

  % Create input data.
  input = ones(N1,N2);

  % Create output data array.
  input = zeros(N1,N2);

  % Broadcast coefficients to everyone else.
  MPI_Bcast( leader, coefs_tag, comm, coefs );

  % Deal input data to everyone else including self.
  for i=1:N2 
    dest = mod((i - 1),comm_size);
    dest_tag = input_tag + i;
    dest_data = input(:,i);
    MPI_Send(dest,dest_tag,comm,dest_data);
  end
end

% Everyone but the leader receives the coefs.
if (my_rank ~= leader)
    % Receive coefs.
    coefs = MPI_Recv( leader, coefs_tag, comm );
end

% Everyone receives the input data and processes the results.
for i=1:N2
  % Compute who the destination is for this message.
  dest = mod((i - 1),comm_size);

  % Check if this destination is me.
  if (my_rank == dest)
    % Compute tags.
    dest_tag = input_tag + i;
    leader_tag = output_tag + i;

    % Receive input.
    i_input =  MPI_Recv(leader,dest_tag,comm);

    % Do computation.
    i_output = ifft(fft(coefs) .* i_input);

    % Send results back to the leader.
    MPI_Send(leader,leader_tag,comm,i_output);
  end
end

% Leader receives all the results.
if (my_rank == leader)
  for i=1:N2 
    % Compute who sent this message.
    dest = mod((i - 1),comm_size);

    leader_tag = output_tag + i;

    % Receive output.
    output(:,i) =  MPI_Recv(dest,leader_tag,comm);
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

