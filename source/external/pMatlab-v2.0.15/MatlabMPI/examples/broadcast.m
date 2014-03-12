%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Matlab MPI script that
% broadcasts a matrix to other processors.
% To run, start Matlab and type:
%
%   eval( MPI_Run('broadcast',4,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('broadcast',4,{'machine1' 'machine2' 'machine3' 'machine4'}) );
%
% Output will be piped into 4 files:
%
%   MatMPI/broadcast.0.out
%   MatMPI/broadcast.1.out
%   MatMPI/broadcast.2.out
%   MatMPI/broadcast.3.out
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

% Set who is source.
source = 0;

% Create a unique tag id for this message (very important in Matlab MPI!).
tag = 1;

% Check to make sure that we more than one processor.
if(comm_size > 1)
  % Source.
  if (my_rank == source)
    % Create data.
    data = 1:10000;
    % Send it.
%    comm = MatMPI_Save_messages(comm,1);
    tic, MPI_Bcast( source, tag, comm, data, data ), toc
%    comm = MatMPI_Save_messages(comm,0);
  end
  % Destination.
  if (my_rank ~= source)
    % Receive data.
    [data data1] = MPI_Recv( source, tag, comm );
    % Check data.
    if(any((data  - (1:10000)) ~= 0))
      disp('ERROR: incorrect data sent.');
      exit;
    end
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

