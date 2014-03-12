%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Matlab MPI script that
% prints out a rank.
%
% To run, start Matlab and type:
%
%   eval( MPI_Run('xbasic',2,{}) );
%
% Or, to run a different machine type:
%
%   eval( MPI_Run('xbasic',2,{'machine1' 'machine2'}) );
%
% Output will be piped into two files:
%
%   MatMPI/xbasic.0.out
%   MatMPI/xbasic.1.out
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

% Modify common directory from default for better performance.
% comm = MatMPI_Comm_dir(comm,'/tmp');

% Get size and rank.
comm_size = MPI_Comm_size(comm);
my_rank = MPI_Comm_rank(comm);

% Print rank.
disp(['my_rank: ',num2str(my_rank)]);

% Wait momentarily.
pause(2.0);

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

